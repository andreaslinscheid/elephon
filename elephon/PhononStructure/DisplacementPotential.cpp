/*	This file DisplacementPotential.cpp is part of elephon.
 *
 *  elephon is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  elephon is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with elephon.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: Jun 30, 2017
 *      Author: A. Linscheid
 */

#include "PhononStructure/DisplacementPotential.h"
#include "LatticeStructure/SymmetryReduction.h"
#include "LatticeStructure/Atom.h"
#include "Algorithms/LinearAlgebraInterface.h"
#include "IOMethods/WriteVASPRealSpaceData.h"
#include "Algorithms/GridRotationMap.h"
#include "Auxillary/UnitConversion.h"
#include <iostream>

namespace elephon
{
namespace PhononStructure
{

void
DisplacementPotential::build(std::shared_ptr<const LatticeStructure::UnitCell> unitCell,
		std::shared_ptr<const LatticeStructure::UnitCell> superCell,
		std::shared_ptr<const std::vector<LatticeStructure::AtomDisplacement>> irredDispl,
		LatticeStructure::RegularBareGrid unitcellGrid,
		LatticeStructure::RegularBareGrid const & supercellGrid,
		std::vector<double> const & potentialUC,
		std::vector< std::vector<double> > potentialDispl,
		std::vector<int> coarseGrainGrid )
{
	unitCell->compute_supercell_dim(superCell, superCellDim_);

	this->set_R_vectors();

	numModes_ = unitCell->get_atoms_list().size()*3;
	const int nrSC = supercellGrid.get_num_points();
	const int nr = nrSC/superCellDim_[0]/superCellDim_[1]/superCellDim_[2];
	nptsRealSpace_ = nr;
	assert( unitcellGrid.get_num_points() == nr );
	assert( potentialDispl.size() == irredDispl->size() );

	// From the potentials, construct the difference
	this->constuct_potential_variation(potentialUC, nrSC, unitcellGrid, supercellGrid, potentialDispl);

	std::vector<std::vector< std::pair<int,std::vector<int> > >> rSuperCellToPrimitve;
	this->build_supercell_to_primitive(
			unitcellGrid,
			supercellGrid,
			unitCell->get_atoms_list(),
			rSuperCellToPrimitve);

	//Regenerate the list of irreducible atoms. This is needed to map the input forces
	// which come for non-equivalent displacements at irreducible atoms into the reducible set.
	std::vector<LatticeStructure::Atom> irredAtoms;
	std::vector<int> redToIrredAtoms, symRedToIrredAtoms;
	std::vector< std::vector<int> > irredToRedAtoms, symIrredToRedAtoms;
	LatticeStructure::SymmetryReduction<LatticeStructure::Atom>(
			unitCell->get_symmetry(),
			unitCell->get_atoms_list(),  irredAtoms,
			redToIrredAtoms, symRedToIrredAtoms,
			irredToRedAtoms, symIrredToRedAtoms);

	//define a consistent atom numbering in the primitive cell.
	//We choose the unitCell.get_atoms_list() order
	std::map< LatticeStructure::Atom, int > primitiveCellLookup;
	for ( int ia = 0 ; ia < unitCell->get_atoms_list().size(); ++ia )
		primitiveCellLookup.insert( std::move( std::make_pair( unitCell->get_atoms_list()[ia], ia ) ) );

	Algorithms::LinearAlgebraInterface linAlg;

	//For each such atom, compute the local displacement potential in x,y and z
	std::vector<double> irreducibleDisplPot( irredAtoms.size()*3*nrSC );

	//This counter advances the point in the set of input potentialVariations.
	//At each irreducible atomic site, we add the irreducible displacements.
	int iIrredDisplCount = 0;
	for ( int ia = 0 ; ia < irredAtoms.size(); ++ia )
	{
		//regenerate displacements to make the connection with the input forces
		LatticeStructure::Symmetry const & sSym = unitCell->get_site_symmetry(ia);

		std::vector<LatticeStructure::AtomDisplacement> irredNotUsed,reducible;
		bool symmetricDispl = not (*irredDispl)[iIrredDisplCount].is_plus_minus_displ_equivalent();
		std::vector<int> redToIrredDispl, symRedToIrredDispl;
		std::vector< std::vector<int> > irredToRedDispl, symIrredToRedDispl;
		unitCell->get_site_displacements( irredAtoms[ia],
				symmetricDispl, sSym,
				1, // not used here
				irredNotUsed,reducible,
				redToIrredDispl, symRedToIrredDispl,
				irredToRedDispl, symIrredToRedDispl);

		//Otherwise the system is underdetermined - this is not supposed to happen.
		int nRedDispl = redToIrredDispl.size();
		assert( nRedDispl >= 3 );

		//build the transpose of the U matrix from the irreducible displacements
		std::vector<double> U( 3 * nRedDispl );
		for ( int idir = 0 ; idir < irredToRedDispl.size();  ++idir)
		{
			auto d = (*irredDispl)[iIrredDisplCount+idir].get_direction();
			for ( auto &xi : d )
				xi *= (*irredDispl)[iIrredDisplCount].get_magnitude();
			for ( int istar = 0 ; istar < symIrredToRedDispl[idir].size(); ++istar)
			{
				auto drot = d;
				sSym.rotate_cartesian( symIrredToRedDispl[idir][istar], drot.begin(), drot.end() );
				std::copy( drot.begin(), drot.end(),  &U[irredToRedDispl[idir][istar]*3] );
			}
		}

		//compute the pseudo inverse
		std::vector<double> piU;
		linAlg.pseudo_inverse( std::move(U), nRedDispl, 3, piU );

		//Obtain the displaced atom in the coordinates of the supercell and compute the rotation maps
		auto p = irredAtoms[ia].get_position();
		for ( int i = 0 ; i < 3 ; ++i)
			p[i] /= double(superCellDim_[i]);
		std::vector< std::vector<int> > rotMapsSiteSymmetry;
		this->compute_rot_map( p, supercellGrid, sSym, rotMapsSiteSymmetry);

		//Symmetry-expand the set of potential variations
		std::vector<double> deltaV( nRedDispl*nrSC );
		for ( int idir = 0 ; idir < irredToRedDispl.size(); ++idir)
		{
			for ( int istar = 0 ; istar < symIrredToRedDispl[idir].size(); ++istar)
			{
				int id = irredToRedDispl[idir][istar];
				assert( id < nRedDispl );
				int isym = symIrredToRedDispl[idir][istar];
				for ( int ir = 0 ; ir < nrSC; ++ir )
				{
					assert( iIrredDisplCount+idir < irredToRedDispl.size() );
					int irRot = rotMapsSiteSymmetry[isym][ir];
					deltaV[id*nrSC+ir] = potentialDispl[iIrredDisplCount+idir][irRot];
				}
			}
		}

		//multiplying the Moore-Penrose pseudo inverse, of the displacements,
		//we obtain the least square fit of transpose of the linear displacement potential
		std::vector<double> linDVscf(3*nrSC);
		linAlg.matrix_matrix_prod( piU, deltaV, linDVscf, 3, nrSC );

		//We convert to a layout with space as the slower running dimension for later
		//symmetry operations and relations. Thus v_[x,y,z](r) is layed out as at
		//point r=r0 we have [v_x(r0),v_y(r0),v_z(r0)]
		for ( int irSC = 0 ; irSC < nrSC; ++irSC)
			for ( int i = 0 ; i < 3; ++i)
				irreducibleDisplPot[(ia*nrSC+irSC)*3+i] = linDVscf[i*nrSC+irSC];

		iIrredDisplCount += irredToRedDispl.size();
	}

	data_.assign(numModes_*nr*this->get_num_R(), 0.0f);
	std::vector<double> linDVscf(3*nrSC);
	std::vector<double> linDVscfTransform(3*nrSC);
	for ( int irA = 0; irA < irredToRedAtoms.size(); ++irA )
	{
		std::copy( &irreducibleDisplPot[irA*3*nrSC],
				   &irreducibleDisplPot[irA*3*nrSC]+3*nrSC,
				   linDVscf.data());

		std::vector<double> gridThisAtom(3*nrSC);
		for ( int ir = 0 ; ir < nrSC; ++ir)
		{
			auto r = supercellGrid.get_vector_direct(ir);
			for ( int xi = 0 ; xi < 3 ; ++ xi)
				gridThisAtom[ir*3+xi] = r[xi];
		}

		//loop the star of the irreducible atom
		for ( int istar = 0; istar < symIrredToRedAtoms[irA].size(); ++istar )
		{
			auto symmetryShiftedGrid = gridThisAtom;
			//transform the linear displacement potential according
			//to the present symmetry operation which takes us from the irreducible
			//to the reducible atom
			int isym = symIrredToRedAtoms[irA][istar];

			//Re-shuffle the field according to the symmetry operation beta = S( alpha )
			//Note: at this point we need that the symmetry of the unit cell is the same as
			//		the one of the supercell
			superCell->get_symmetry().apply(isym,symmetryShiftedGrid.begin(),symmetryShiftedGrid.end(),true);
			//since the grid is regular, we must obtain the (x,y,z) indices by scaling with the grid dimension
			//in each direction
			std::vector<int> xyz(3);
			for ( int ir = 0 ; ir < nrSC; ++ir)
			{
				for ( int j = 0 ; j < 3; ++j)
				{
					symmetryShiftedGrid[ir*3+j] -= std::floor(symmetryShiftedGrid[ir*3+j]);
					symmetryShiftedGrid[ir*3+j] *= supercellGrid.get_grid_dim()[j];
					xyz[j] = std::floor(symmetryShiftedGrid[ir*3+j]+0.5);
					if ( std::abs(symmetryShiftedGrid[ir*3+j]-xyz[j]) > 0.01 )
						throw std::logic_error("The symmetry operation that transports the atom "
								"does not map the grid to itself which can't be.");
				}
				int cnsq = supercellGrid.get_xyz_to_reducible(xyz);
				assert( (cnsq >= 0) && (cnsq < nrSC) );

				// re-shuffel the vector field to the right location in the grid without transforming it yet.
				for (int xi = 0 ; xi < 3 ; ++xi)
					linDVscfTransform[cnsq*3+xi]=linDVscf[ir*3+xi];
			}
			// transform the vector field.
			superCell->get_symmetry().rotate_cartesian(
					isym,
					linDVscfTransform.begin(),
					linDVscfTransform.end());

			//Convert into the format of the linear derivative potential in real space
			int iaDispl = irredToRedAtoms[irA][istar];
			for ( int ir = 0 ; ir < nrSC; ++ir)
			{
				for ( int xi_displ = 0 ; xi_displ < 3; ++xi_displ)
				{
					int mu = iaDispl*3 + xi_displ;
					int ir1uc = rSuperCellToPrimitve[iaDispl][ir].first;
					std::vector<int> R = rSuperCellToPrimitve[iaDispl][ir].second;
					for ( int i = 0 ; i < 3 ; ++i)
						R[i] = R[i] < 0 ? R[i] + RVectorDim_[i] : R[i];
					int iR = this->RVectorLayout(R[0],R[1],R[2]);
					data_[this->mem_layout(ir1uc,mu,iR)] = linDVscfTransform[ir*3+xi_displ];
				}
			}
		}
	}

	//keep a copy for lattice information, symmetry ect ...
	unitCell_ = unitCell;
	unitCellGrid_ = std::move(unitcellGrid);

	this->clean_displacement_potential();

	if (not coarseGrainGrid.empty())
	{
		assert(coarseGrainGrid.size() == 3);
		coarseGrainGrid_ = std::make_shared<LatticeStructure::RegularBareGrid>();
		coarseGrainGrid_->initialize(
				std::move(coarseGrainGrid),
				/*is reciprcal space=*/true,
				unitCellGrid_.get_grid_prec(),
				/*set grid shit to zero :*/{0.0, 0.0, 0.0},
				unitcellGrid.get_lattice());
	}
}

void
DisplacementPotential::compute_rot_map(
		std::vector<double> const & shift,
		LatticeStructure::RegularBareGrid const & supercellGrid,
		LatticeStructure::Symmetry const & siteSymmetry,
		std::vector< std::vector<int> > & rotMap) const
{
	Algorithms::compute_grid_rotation_map(shift, supercellGrid, siteSymmetry, rotMap);
}

void
DisplacementPotential::compute_dvscf_q(
		std::vector<double> const & qVect,
		std::vector<double> const & modes,
		std::vector<std::complex<double>> const & dynamicalMatrices,
		std::vector<double> const & masses,
		std::vector<std::complex<float>> & dvscf,
		std::vector<std::vector<std::complex<float>>> & buffer,
		double freqCutoff) const
{
	assert( qVect.size() % 3 == 0 );
	const int nq = qVect.size()/3;
	const int nR = this->get_num_R();
	const int nAUC = numModes_/3;
	assert( modes.size() == numModes_*nq );
	assert( (dynamicalMatrices.size()/nq) / (numModes_*numModes_) == 1 );
	assert( (masses.size()*3) / numModes_ == 1 );
	freqCutoff = freqCutoff < 1e-5 ? 1e-5 : freqCutoff;

	// introduce names for the various buffers
	buffer.resize(2);
	std::vector<std::complex<float>> & ftDisplPot = buffer[0];
	ftDisplPot.resize(numModes_*nptsRealSpace_);
	std::vector<std::complex<float>> & dynmatMassBuffer = buffer[1];

	dvscf.assign(nq*numModes_*nptsRealSpace_, 0.0f);
	Algorithms::LinearAlgebraInterface linAlg;
	for ( int iq = 0 ; iq < nq; ++iq)
	{
		std::fill(ftDisplPot.begin(), ftDisplPot.end(), std::complex<float>(0));
		for (int ia = 0 ; ia < nAUC; ++ia)
		{
			double tau[3] = {unitCell_->get_atoms_list()[ia].get_position()[0],
							 unitCell_->get_atoms_list()[ia].get_position()[1],
							 unitCell_->get_atoms_list()[ia].get_position()[2] };
			for ( int iR = 0 ; iR < nR; ++iR)
			{
				float dprod = -2.0*M_PI*(qVect[iq*3+0]*(RVectors_[3*iR+0] + tau[0])
										+qVect[iq*3+1]*(RVectors_[3*iR+1] + tau[1])
										+qVect[iq*3+2]*(RVectors_[3*iR+2] + tau[2]));
				std::complex<float> phase = std::exp( std::complex<float>(0,dprod));
				for (int iDispl = 0 ; iDispl < 3; ++iDispl)
				{
					int mu = ia*3+iDispl;
					for (int ir = 0 ; ir < nptsRealSpace_; ++ir)
						ftDisplPot[mu*nptsRealSpace_+ir] += phase*data_[this->mem_layout(ir, mu, iR)];
				}
			}
		}

		// in case there are only few modes in the system, avoid the calling overhead of
		// the blas package
		if (nAUC < 3)
		{
			for ( int mu1 = 0 ; mu1 < numModes_; ++mu1 )
				for ( int mu2 = 0 ; mu2 < numModes_; ++mu2 )
					for ( int ir = 0 ; ir < nptsRealSpace_; ++ir )
						dvscf[(iq*numModes_+mu1)*nptsRealSpace_+ir] +=
								std::complex<float>(dynamicalMatrices[(iq*numModes_+mu1)*numModes_+mu2]) / float(masses[mu2/3])
										  *ftDisplPot[mu2*nptsRealSpace_+ir];
		}
		else
		{
			dynmatMassBuffer.resize(numModes_*numModes_);
			for ( int mu1 = 0 ; mu1 < numModes_; ++mu1 )
				for ( int mu2 = 0 ; mu2 < numModes_; ++mu2 )
					dynmatMassBuffer[mu1*numModes_+mu2] =
							std::complex<float>(dynamicalMatrices[(iq*numModes_+mu1)*numModes_+mu2]) / float(masses[mu2/3]);

			auto r_ptr = &dvscf[iq*numModes_*nptsRealSpace_];
			linAlg.call_gemm( 'n', 'n',
					numModes_, nptsRealSpace_ , numModes_,
					std::complex<float>(1.0f), dynmatMassBuffer.data(), numModes_,
					ftDisplPot.data(), nptsRealSpace_,
					std::complex<float>(0.0f), r_ptr, nptsRealSpace_);
		}

		for ( int mu = 0 ; mu < numModes_; ++mu )
		{
			float scale = (std::abs(modes[iq*numModes_+mu]) < freqCutoff ? 0.0f :
					static_cast<float>(Auxillary::units::SQRT_HBAR_BY_2M_THZ_TO_ANGSTROEM / std::sqrt(modes[iq*numModes_+mu])) );
			for (int ir = 0 ; ir < nptsRealSpace_; ++ir)
				dvscf[(iq*numModes_+mu)*nptsRealSpace_+ir] *= scale;
		}
	}
}

void
DisplacementPotential::query_q(
		std::vector<double> const & qVect,
		std::map<LatticeStructure::RegularBareGrid::GridPoint, std::vector<int>> & qGridCorseGainToQIndex) const
{
	assert(qVect.size()%3 == 0);
	qGridCorseGainToQIndex.clear();
	const int nq = qVect.size()/3;
	const double gPrec = unitCell_->get_symmetry().get_symmetry_prec();

	if ( not coarseGrainGrid_ )
	{
		for (int iq = 0 ; iq < nq ;++iq)
		{
			LatticeStructure::RegularBareGrid::GridPoint g(std::move(std::vector<double>(&qVect[iq*3], &qVect[iq*3]+3)), gPrec);
			auto ret = qGridCorseGainToQIndex.insert(std::move(std::make_pair(std::move(g), std::vector<int>())));
			ret.first->second.push_back(iq);
		}
	}
	else
	{
		std::vector<int> closestGridIndices;
		coarseGrainGrid_->find_closest_reducible_grid_points(qVect, closestGridIndices);

		for (int iq = 0 ; iq < nq ;++iq)
		{
			auto gridVector = coarseGrainGrid_->get_vector_direct(closestGridIndices[iq]);
			LatticeStructure::RegularBareGrid::GridPoint g(std::move(gridVector), gPrec);
			auto ret = qGridCorseGainToQIndex.insert(std::move(std::make_pair(std::move(g), std::vector<int>())));
			ret.first->second.push_back(iq);
		}
	}
}

int
DisplacementPotential::get_num_R() const
{
	return RVectorDim_[0]*RVectorDim_[1]*RVectorDim_[2];
}

int
DisplacementPotential::RVectorLayout(int iRx, int iRy, int iRz) const
{
	assert((iRx >= 0) && (iRx<RVectorDim_[0]));
	assert((iRy >= 0) && (iRy<RVectorDim_[1]));
	assert((iRz >= 0) && (iRz<RVectorDim_[2]));
	return (iRz*RVectorDim_[1]+iRy)*RVectorDim_[0]+iRx;
}

int
DisplacementPotential::get_num_modes() const
{
	return numModes_;
}

int
DisplacementPotential::mem_layout(int ir, int mu, int iR ) const
{
	assert(ir < this->nptsRealSpace_ );
	assert(mu < this->numModes_);
	assert(iR < this->get_num_R());
	return ir+nptsRealSpace_*(mu+numModes_*iR);
}

LatticeStructure::RegularBareGrid const &
DisplacementPotential::get_real_space_grid() const
{
	return unitCellGrid_;
}

void
DisplacementPotential::write_dvscf(
		int atomIndex, int xi,
		std::string filename) const
{
	assert((atomIndex >= 0) && (atomIndex<numModes_/3));
	assert((xi >= 0) && (xi < 3));
	IOMethods::WriteVASPRealSpaceData writer;
	std::string comment = "Displacement potential in real space; atom # "+ std::to_string(atomIndex) +", vibration in dir. "
			+ (xi == 0 ? "x" : (xi == 1 ? "y" : "z") ) + "\n";

	LatticeStructure::UnitCell supercell;

	// build a supercell where the displaced atom is in the center.
	std::vector<LatticeStructure::Atom> scAtoms;
	scAtoms.reserve(unitCell_->get_atoms_list().size()*superCellDim_[0]*superCellDim_[1]*superCellDim_[2]);
	auto atomsUC = unitCell_->get_atoms_list();
	auto tau0 = atomsUC[atomIndex].get_position();
	for (auto i = 0 ; i < tau0.size(); ++i)
		tau0[i] /=  superCellDim_[i];

	for (int isx = 0 ; isx < superCellDim_[0]; ++isx)
		for (int isy = 0 ; isy < superCellDim_[1]; ++isy)
			for (int isz = 0 ; isz < superCellDim_[2]; ++isz)
			{
				int is[] = {isx, isy, isz};
				for (auto a: atomsUC)
				{
					auto pos = a.get_position();
					for (auto i = 0 ; i < pos.size(); ++i)
						pos[i] = (is[i]+pos[i])/superCellDim_[i] - tau0[i] + 0.5; // We shift the entire system by 0.5 to have is displayed from 0-1
					a.set_position(std::move(pos));
					scAtoms.push_back(std::move(a));
				}
			}
	supercell.initialize(	std::move(scAtoms),
							unitCell_->get_lattice().build_supercell(superCellDim_),
							unitCell_->get_symmetry());

	auto dim = unitCellGrid_.get_grid_dim();
	for (int i = 0 ; i < 3 ; ++i)
		dim[i] *= superCellDim_[i];

	LatticeStructure::RegularBareGrid scGrid;
	scGrid.initialize(dim,
			false,
			unitCellGrid_.get_grid_prec(),
			{0.0, 0.0, 0.0},
			supercell.get_lattice() );

	std::vector<std::vector< std::pair<int,std::vector<int> > >> rSuperCellToPrimitve;
	this->build_supercell_to_primitive( unitCellGrid_, scGrid, unitCell_->get_atoms_list(), rSuperCellToPrimitve);

	std::vector<double> dataSC( scGrid.get_num_points() );
	assert( rSuperCellToPrimitve[atomIndex].size() == scGrid.get_num_points() );
	for ( int irSC = 0 ; irSC < scGrid.get_num_points(); ++irSC )
	{
		auto ir = rSuperCellToPrimitve[atomIndex][irSC].first;
		auto & R = rSuperCellToPrimitve[atomIndex][irSC].second;
		int iRx = R[0] < 0 ? R[0] + RVectorDim_[0] : R[0];
		int iRy = R[1] < 0 ? R[1] + RVectorDim_[1] : R[1];
		int iRz = R[2] < 0 ? R[2] + RVectorDim_[2] : R[2];
		int iR = this->RVectorLayout(iRx, iRy, iRz);

		// shift the data position by 0.5 in each direction
		auto xyz = scGrid.get_reducible_to_xyz(irSC);
		for (int i = 0 ; i < xyz.size() ; ++i)
			xyz[i] -= scGrid.get_grid_dim()[i]/2;
		int irSC_shift = scGrid.get_xyz_to_reducible_periodic(xyz);

		dataSC[irSC_shift] = data_[this->mem_layout( ir , atomIndex*3+xi, iR)];
	}

	writer.write_file(	filename,
						comment,
						scGrid.get_grid_dim(),
						std::make_shared<decltype(supercell)>(std::move(supercell)),
						dataSC );
}

void
DisplacementPotential::write_dvscf_q(
		std::vector<double> const & qVect,
		std::vector<int> modeIndices,
		std::vector<double> const & modes,
		std::vector<std::complex<double>> const & dynamicalMatrices,
		std::vector<double> const & masses,
		std::string filename) const
{
	assert( qVect.size()%3 == 0 );
	int nq = qVect.size()/3;

	std::vector<std::complex<float>> qDataPlus;
	std::vector<std::vector<std::complex<float>>> buffers;
	this->compute_dvscf_q(qVect, modes, dynamicalMatrices, masses, qDataPlus, buffers);

	int nr = unitCellGrid_.get_num_points();
	assert(qDataPlus.size() == nq*numModes_*nr);

	auto mod_filename = [] (std::string filename, int iq, int mu) {
		if ( filename.length() >= 5 )
			if ( filename.substr(filename.length()-4,filename.length()).compare( ".dat" ) == 0 )
				filename.insert(filename.length()-4,"_"+std::to_string(iq)+"_"+std::to_string(mu));
		return filename;
	};

	IOMethods::WriteVASPRealSpaceData writer;
	for ( int iq = 0; iq < nq ; ++iq )
	{
		for ( auto mu : modeIndices )
		{
			assert( (mu >= 0) && (mu < numModes_) );
			std::string comment = "Displacement potential dvscf(r)/du(q,mu); q = ("
					+std::to_string(qVect[iq*3+0])+" , "+std::to_string(qVect[iq*3+1])+" , "+std::to_string(qVect[iq*3+2])
					+") mode # "+ std::to_string(mu) +". Units are eV.\n";
			std::vector<double> data(nr);
			for ( int ir = 0 ; ir < nr ; ++ir )
				data[ir] = std::real(qDataPlus[(iq*numModes_+mu)*nr+ir]);
			writer.write_file( mod_filename(filename,iq,mu), comment,
					unitCellGrid_.get_grid_dim(), unitCell_, data );
		}
	}
}

void
DisplacementPotential::build_supercell_to_primitive(
		LatticeStructure::RegularBareGrid const & primitiveCellGrid,
		LatticeStructure::RegularBareGrid const & supercellGrid,
		std::vector<LatticeStructure::Atom> const & atomsUC,
		std::vector< std::vector< std::pair<int,std::vector<int> > > > & rSuperCellToPrimitve) const
{
	const int nRSC = supercellGrid.get_num_points();

	std::vector<double> rVectorsSupercell(nRSC*3);
	for ( int irSC = 0 ; irSC < nRSC ; ++irSC )
	{
		auto rvec = supercellGrid.get_vector_direct(irSC);
		assert(rvec.size() == 3);
		for ( int i = 0 ; i < 3 ; ++i )
			rVectorsSupercell[irSC*3+i] = rvec[i];
	}

	rSuperCellToPrimitve.resize(atomsUC.size());
	for ( int ia = 0 ; ia < atomsUC.size() ; ++ia )
	{
		auto aPos = atomsUC[ia].get_position();
		rSuperCellToPrimitve[ia].resize(nRSC);

		for ( int irSC = 0 ; irSC < nRSC ; ++irSC )
		{
			std::vector<int> R(3);
			std::vector<int> xyz(3);
			for ( int i = 0 ; i < 3 ; ++i )
			{
				double rvec = rVectorsSupercell[irSC*3+i];

				// scale to units of the primitive cell
				// todo use map_supercell_to_primitive();
				rvec *= superCellDim_[i];

				// transform to a coordinate system, where the current atom is in the center
				// NOTE: the atom coordinates are not necessarily on the grid.
				rvec -= aPos[i];

				// obtain the vector in the vector in the first unit cell and the Lattice vector taking it to the position
				// in the supercell
				R[i] = std::floor(rvec+0.5);
				double rvecUC = rvec-R[i];
				assert( (rvecUC >= -0.5-unitCellGrid_.get_grid_prec()) && (rvecUC < 0.5+unitCellGrid_.get_grid_prec()));
				if ( rvecUC < -0.5)
					rvecUC = -0.5;
				if ( rvecUC > 0.5 )
					rvecUC = 0.5;

				//Convert to a coordinate index in the range [0,1[ and fetch the grid index
				rvecUC -= std::floor(rvecUC);
				rvecUC *= primitiveCellGrid.get_grid_dim()[i];
				xyz[i] = std::floor(rvecUC+0.5);
			}
			// NOTE: since the atom position shit means we are out of the grid, we have to allow for
			// the case where any xyz[i] is equal to primitiveCellGrid.get_grid_dim()[i]. Thus we have to use the
			// periodic version get_xyz_to_reducible.
			int cnsq = primitiveCellGrid.get_xyz_to_reducible_periodic(xyz);
			assert( (cnsq >= 0) and ( cnsq < primitiveCellGrid.get_num_points() ) );
			rSuperCellToPrimitve[ia][irSC] = std::move(std::make_pair(cnsq, std::move(R) ) );
		}
	}
}

void
DisplacementPotential::clean_displacement_potential()
{
	const int numAtoms = this->get_num_modes()/3;
	const int nr = unitCellGrid_.get_num_points();
	const int nR = this->get_num_R();

	// apply the sum rule
	for ( int xi = 0 ; xi < 3 ; ++xi)
	{
		double integral = 0.0;
		for (int iR = 0; iR < nR ; ++iR)
			for (int ia = 0 ; ia < numAtoms; ++ia)
				for (int ir = 0; ir < nr ; ++ir)
					integral += data_[this->mem_layout( ir , ia*3+xi, iR)];
		integral /= static_cast<double>(nr*nR*numAtoms);

		for (int iR = 0; iR < nR ; ++iR)
			for (int ia = 0 ; ia < numAtoms; ++ia)
				for (int ir = 0; ir < nr ; ++ir)
					data_[this->mem_layout( ir , ia*3+xi, iR)] -= integral;
	}
}

void
DisplacementPotential::set_R_vectors()
{
	RVectorDim_.resize(3);
	for (int i = 0 ; i < 3; ++i)
		RVectorDim_[i] = superCellDim_[i] % 2 == 0 ? superCellDim_[i]+1 : superCellDim_[i];

	RVectors_.resize(3*RVectorDim_[0]*RVectorDim_[1]*RVectorDim_[2]);
	for ( int iRz = 0 ; iRz < RVectorDim_[2]; ++iRz )
		for ( int iRy = 0 ; iRy < RVectorDim_[1]; ++iRy )
			for ( int iRx = 0 ; iRx < RVectorDim_[0]; ++iRx )
			{
				int iR = this->RVectorLayout(iRx, iRy, iRz);
				int R[] = {	iRx <= RVectorDim_[0]/2 ? iRx : iRx - RVectorDim_[0],
							iRy <= RVectorDim_[1]/2 ? iRy : iRy - RVectorDim_[1],
							iRz <= RVectorDim_[2]/2 ? iRz : iRz - RVectorDim_[2] };
				std::copy(R, R+3, &RVectors_[iR*3]);
			}
}

void
DisplacementPotential::constuct_potential_variation(
		std::vector<double> const & potentialUC,
		int nRSC,
		LatticeStructure::RegularBareGrid const & unitcellGrid,
		LatticeStructure::RegularBareGrid const & supercellGrid,
		std::vector< std::vector<double> > & potentialDispl)
{
	std::vector<int> mapback(nRSC);
	std::vector<double> rVectorsSupercell(nRSC*3);
	for ( int irSC = 0 ; irSC < nRSC ; ++irSC )
	{
		std::vector<int> xyz = supercellGrid.get_reducible_to_xyz(irSC);
		for (int i = 0 ; i < 3 ; ++i)
			xyz[i] = xyz[i]%unitcellGrid.get_grid_dim()[i];
		int irUC = unitcellGrid.get_xyz_to_reducible(xyz);
		assert((irUC >= 0) && (irUC < unitcellGrid.get_num_points()));
		mapback[irSC] = irUC;
	}

	for (int ird = 0 ; ird < potentialDispl.size(); ++ird)
	{
		assert(potentialDispl[ird].size() == nRSC);
		for (int ir = 0 ; ir < nRSC; ++ir)
			potentialDispl[ird][ir] -= potentialUC[mapback[ir]];
	}
}

} /* namespace PhononStructure */
} /* namespace elephon */
