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
		std::vector< std::vector<double> > const & potentialDispl,
		std::vector<int> coarseGrainGrid )
{
	unitCell->compute_supercell_dim(superCell, superCellDim_);
	numModes_ = unitCell->get_atoms_list().size()*3;
	int nRSC = supercellGrid.get_num_points();
	int nR = nRSC/this->get_num_R();
	nptsRealSpace_ = nR;
	assert( unitcellGrid.get_num_points() == nR );
	assert( potentialDispl.size() == irredDispl->size() );

	std::vector< std::pair<int,std::vector<int> > > rSuperCellToPrimitve;
	this->build_supercell_to_primite(unitcellGrid, supercellGrid, rSuperCellToPrimitve);

	//From the potentials, construct the difference
	std::vector< std::vector<double> > potentialVariation = potentialDispl;
	for (int ird = 0 ; ird < potentialDispl.size(); ++ird)
	{
		assert(potentialVariation[ird].size() == nRSC);
		for (int ir = 0 ; ir < nRSC; ++ir)
		{
			int firstUC = rSuperCellToPrimitve[ir].first;
			assert( (firstUC >= 0) and (firstUC < nR) );
			potentialVariation[ird][ir] -= potentialUC[firstUC];
		}
	}

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
	std::vector<double> irreducibleDisplPot( irredAtoms.size()*3*nRSC );

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
		std::vector<double> deltaV( nRedDispl*nRSC );
		for ( int idir = 0 ; idir < irredToRedDispl.size(); ++idir)
		{
			for ( int istar = 0 ; istar < symIrredToRedDispl[idir].size(); ++istar)
			{
				int id = irredToRedDispl[idir][istar];
				assert( id < nRedDispl );
				int isym = symIrredToRedDispl[idir][istar];
				for ( int ir = 0 ; ir < nRSC; ++ir )
				{
					assert( iIrredDisplCount+idir < irredToRedDispl.size() );
					int irRot = rotMapsSiteSymmetry[isym][ir];
					deltaV[id*nRSC+ir] = potentialVariation[iIrredDisplCount+idir][irRot];
				}
			}
		}

		//multiplying the Moore-Penrose pseudo inverse, of the displacements,
		//we obtain the least square fit of transpose of the linear displacement potential
		std::vector<double> linDVscf(3*nRSC);
		linAlg.matrix_matrix_prod( piU, deltaV, linDVscf, 3, nRSC );

		//We convert to a layout with space as the slower running dimension for later
		//symmetry operations and relations. Thus v_[x,y,z](r) is layed out as at
		//point r=r0 we have [v_x(r0),v_y(r0),v_z(r0)]
		for ( int iaSC = 0 ; iaSC < nRSC; ++iaSC)
			for ( int i = 0 ; i < 3; ++i)
				irreducibleDisplPot[(ia*nRSC+iaSC)*3+i] = linDVscf[i*nRSC+iaSC];

		iIrredDisplCount += irredToRedDispl.size();
	}

	data_.resize( this->get_num_R()*numModes_*nR );
	std::vector<double> linDVscf(3*nRSC);
	std::vector<double> linDVscfTransform(3*nRSC);
	for ( int irA = 0; irA < irredToRedAtoms.size(); ++irA )
	{
		std::copy( &irreducibleDisplPot[irA*3*nRSC],
				   &irreducibleDisplPot[irA*3*nRSC]+3*nRSC,
				   linDVscf.data());

		std::vector<double> gridThisAtom(3*nRSC);
		for ( int ir = 0 ; ir < nRSC; ++ir)
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
			for ( int ir = 0 ; ir < nRSC; ++ir)
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
				assert( (cnsq >= 0) && (cnsq < nRSC) );

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
			for ( int ir = 0 ; ir < nRSC; ++ir)
			{
				for ( int xi_displ = 0 ; xi_displ < 3; ++xi_displ)
				{
					int mu = iaDispl*3 + xi_displ;
					int ir1uc = rSuperCellToPrimitve[ir].first;
					std::vector<int> R = rSuperCellToPrimitve[ir].second;
					for ( int i = 0 ; i < 3 ; ++i)
						R[i] = R[i] < 0 ? R[i] + superCellDim_[i] : R[i];
					int iR = this->RVectorLayout(R[2],R[1],R[0]);
					data_[this->mem_layout(ir1uc,mu,iR)] = linDVscfTransform[ir*3+xi_displ];
				}
			}
		}
	}

	//keep a copy for lattice information, symmetry ect ...
	unitCell_ = unitCell;
	unitCellGrid_ = std::move(unitcellGrid);

	this->clean_displacement_potential();

	if (RVectors_.size() != nR)
	{
		RVectors_.resize(3*nR);
		for ( int iRz = 0 ; iRz < superCellDim_[2]; ++iRz )
			for ( int iRy = 0 ; iRy < superCellDim_[1]; ++iRy )
				for ( int iRx = 0 ; iRx < superCellDim_[0]; ++iRx )
				{
					int iR = this->RVectorLayout(iRx,iRy,iRz);
					int R[] = {	iRx <= superCellDim_[0]/2 ? iRx : iRx - superCellDim_[0],
									iRy <= superCellDim_[1]/2 ? iRy : iRy - superCellDim_[1],
									iRz <= superCellDim_[2]/2 ? iRz : iRz - superCellDim_[2] };
					std::copy(R, R+3, &RVectors_[iR*3]);
				}
	}

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
		std::vector<std::complex<double>> const & dynamicalMatrices,
		std::vector<double> const & masses,
		std::vector<std::complex<float>> & dvscf,
		std::vector<std::vector<std::complex<float>>> & buffer) const
{
	assert( qVect.size() % 3 == 0 );
	const int nq = qVect.size()/3;
	const int nR = this->get_num_R();
	assert( (dynamicalMatrices.size()/nq) / (numModes_*numModes_) == 1 );
	assert( (masses.size()*3) / numModes_ == 1 );

	// introduce names for the various buffers
	buffer.resize(4);
	std::vector<std::complex<float>> & ftDisplPot = buffer[0];
	std::vector<std::complex<float>> & dvscfqBuffer = buffer[1];
	std::vector<std::complex<float>> & dynmatMassBuffer = buffer[2];
	std::vector<std::complex<float>> & phaseBuffer = buffer[3];

	// in case there are only few modes in the system, avoid the caling overhead of
	// the blas package
	if ( numModes_ < 12 )
	{
		this->compute_dvscf_q_optimized_low_num_modes( qVect, dynamicalMatrices, masses, dvscf, ftDisplPot);
		return;
	}

	Algorithms::LinearAlgebraInterface linAlg;
	// fill the buffers
	if ( dvscfqBuffer.empty() )
		dvscfqBuffer.assign(data_.begin(), data_.end());

	if ( ftDisplPot.size() != nq*numModes_*nptsRealSpace_ )
		ftDisplPot.resize(nq*numModes_*nptsRealSpace_);

	if ( dynmatMassBuffer.size() != numModes_*numModes_)
		dynmatMassBuffer.resize(numModes_*numModes_);

	phaseBuffer.resize(nq*nR);
	for ( int iq = 0 ; iq < nq; ++iq)
	{
		for ( int iR = 0 ; iR < nR; ++iR)
		{
			float dprod = 2.0*M_PI*(qVect[iq*3+0]*RVectors_[3*iR+0]
									+qVect[iq*3+1]*RVectors_[3*iR+1]
									+qVect[iq*3+2]*RVectors_[3*iR+2]);
			phaseBuffer[iq*nR+iR] = std::exp( std::complex<float>(0,dprod));
		}
	}

	dvscf.resize(nq*numModes_*nptsRealSpace_);

	// this formulates the Fourier transform as a matrix multiplication (sum convention)
	// dvscf(iq,[ir,imu]) = phase(iq,iR)*dvscf(iR,[ir,imu])
	// where the brackets are in column major format.

	linAlg.call_gemm( 'n', 'n',
			nq, nptsRealSpace_*numModes_, nR,
			std::complex<float>(1.0f), phaseBuffer.data(), nR,
			dvscfqBuffer.data(), numModes_*nptsRealSpace_,
			std::complex<float>(0.0f), ftDisplPot.data(), nptsRealSpace_*numModes_);

	for ( int iq = 0 ; iq < nq; ++iq)
	{
		for ( int mu1 = 0 ; mu1 < numModes_; ++mu1 )
			for ( int mu2 = 0 ; mu2 < numModes_; ++mu2 )
				dynmatMassBuffer[mu1*numModes_+mu2] = dynamicalMatrices[(iq*numModes_+mu1)*numModes_+mu2] / masses[mu2/3];

		auto r_ptr = &dvscf[iq*numModes_*nptsRealSpace_];
		auto dmat_ptr = &dynmatMassBuffer[0];
		auto ftdp_ptr = &ftDisplPot[iq*numModes_*nptsRealSpace_];
		linAlg.call_gemm( 'n', 'n',
				numModes_, nptsRealSpace_ , numModes_,
				std::complex<float>(1.0f), dmat_ptr, numModes_,
				ftdp_ptr, nptsRealSpace_,
				std::complex<float>(0.0f), r_ptr, nptsRealSpace_);
	}
}

void
DisplacementPotential::compute_dvscf_q_optimized_low_num_modes(
		std::vector<double> const & qVect,
		std::vector<std::complex<double>> const & dynamicalMatrices,
		std::vector<double> const & masses,
		std::vector<std::complex<float>> & dvscf,
		std::vector<std::complex<float>> & ftDisplPot) const
{
	const int nq = qVect.size()/3;
	const int nR = this->get_num_R();

	if ( ftDisplPot.size() != numModes_*nptsRealSpace_ )
		ftDisplPot.resize(numModes_*nptsRealSpace_);

	dvscf.assign(nq*numModes_*nptsRealSpace_, std::complex<float>(0));
	for ( int iq = 0 ; iq < nq; ++iq)
	{
		std::fill(ftDisplPot.begin(), ftDisplPot.end(), std::complex<float>(0));
		for ( int iR = 0 ; iR < nR; ++iR)
		{
			float dprod = 2.0*M_PI*(qVect[iq*3+0]*RVectors_[3*iR+0]
												+qVect[iq*3+1]*RVectors_[3*iR+1]
												+qVect[iq*3+2]*RVectors_[3*iR+2]);
			std::complex<float> phase = std::exp( std::complex<float>(0,dprod));

			for (int imr = 0 ; imr < nptsRealSpace_*numModes_; ++imr)
				ftDisplPot[imr] += phase*data_[imr+nptsRealSpace_*numModes_*iR];
		}

		for ( int mu1 = 0 ; mu1 < numModes_; ++mu1 )
			for ( int mu2 = 0 ; mu2 < numModes_; ++mu2 )
				for ( int ir = 0 ; ir < nptsRealSpace_; ++ir )
					dvscf[(iq*numModes_+mu1)*nptsRealSpace_+ir] +=
							std::complex<float>(dynamicalMatrices[(iq*numModes_+mu1)*numModes_+mu2]) / float(masses[mu2/3])
									  *ftDisplPot[mu2*nptsRealSpace_+ir];
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
	return superCellDim_[0]*superCellDim_[1]*superCellDim_[2];
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

int
DisplacementPotential::RVectorLayout(int iRz, int iRy, int iRx ) const
{
	assert( (iRx >= 0) && (iRx < superCellDim_[0]) );
	assert( (iRy >= 0) && (iRy < superCellDim_[1]) );
	assert( (iRz >= 0) && (iRz < superCellDim_[2]) );
	return (iRz*superCellDim_[1]+iRy)*superCellDim_[0]+iRx;
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
	//This method write the potential in the supercell
	auto sc = unitCell_->build_supercell( superCellDim_[0], superCellDim_[1], superCellDim_[2] );
	elephon::LatticeStructure::RegularBareGrid scGrid;
	auto dim = unitCellGrid_.get_grid_dim();
	scGrid.initialize(
			std::vector<int>({dim[0]*superCellDim_[0],dim[1]*superCellDim_[1],dim[2]*superCellDim_[2]}),
			false,
			unitCellGrid_.get_grid_prec(),
			unitCellGrid_.get_grid_shift(),
			sc.get_lattice() );

	std::vector< std::pair<int,std::vector<int> > > rSuperCellToPrimitve;
	this->build_supercell_to_primite( unitCellGrid_, scGrid, rSuperCellToPrimitve);

	std::vector<double> dataSC( scGrid.get_num_points() );
	assert( rSuperCellToPrimitve.size() == scGrid.get_num_points() );
	for ( int irSC = 0 ; irSC < scGrid.get_num_points(); ++irSC )
	{
		auto ir = rSuperCellToPrimitve[irSC].first;
		auto & R = rSuperCellToPrimitve[irSC].second;
		int iRx = R[0] < 0 ? R[0] + superCellDim_[0] : R[0];
		int iRy = R[1] < 0 ? R[1] + superCellDim_[1] : R[1];
		int iRz = R[2] < 0 ? R[2] + superCellDim_[2] : R[2];
		dataSC[irSC] = data_[this->mem_layout( ir , atomIndex*3+xi, this->RVectorLayout(iRz,iRy,iRx))];
	}
	writer.write_file(filename, comment, scGrid.get_grid_dim(), std::make_shared<decltype(sc)>(std::move(sc)), dataSC );
}

void
DisplacementPotential::write_dvscf_q(
		std::vector<double> const & qVect,
		std::vector<int> modeIndices,
		std::vector<std::complex<double>> const & dynamicalMatrices,
		std::vector<double> const & masses,
		std::string filename) const
{
	assert( qVect.size()%3 == 0 );
	int nq = qVect.size()/3;

	std::vector<std::complex<float>> qDataPlus;
	std::vector<std::vector<std::complex<float>>> buffers;
	this->compute_dvscf_q(qVect, dynamicalMatrices, masses, qDataPlus, buffers);

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
					+") mode # "+ std::to_string(mu) +"\n";
			std::vector<double> data(nr);
			for ( int ir = 0 ; ir < nr ; ++ir )
				data[ir] = std::real(qDataPlus[(iq*numModes_+mu)*nr+ir]);
			writer.write_file( mod_filename(filename,iq,mu), comment,
					unitCellGrid_.get_grid_dim(), unitCell_, data );
		}
	}
}

void
DisplacementPotential::build_supercell_to_primite(
		LatticeStructure::RegularBareGrid const & primitiveCellGrid,
		LatticeStructure::RegularBareGrid const & supercellGrid,
		std::vector< std::pair<int,std::vector<int> > > & rSuperCellToPrimitve) const
{
	const int nRSC = supercellGrid.get_num_points();

	//Create a table which maps a point in real space in the supercell
	// into a point in the primitive cell plus a lattice vector
	rSuperCellToPrimitve.resize( nRSC );
	for ( int irSC = 0 ; irSC < nRSC ; ++irSC )
	{
		auto rvec = supercellGrid.get_vector_direct(irSC);
		std::vector<int> R(3);
		std::vector<double> rvecUC(3);
		std::vector<int> xyz(3);
		for ( int i = 0 ; i < 3 ; ++i )
		{
			const double g = primitiveCellGrid.get_grid_prec();
			//scale to units of the primitive cell
			rvec[i] *= superCellDim_[i];
			//We add a tiny bit to the vector so that 1.9999999 will not (incorrectly) map
			//to R = 1 and then 0.9999999 will (correctly) map to the lattice site at 0.
			R[i] = std::floor(rvec[i]+0.5+2*g);
			rvecUC[i] = rvec[i]-R[i];
			assert( (rvecUC[i] >= -0.5-g) && (rvecUC[i] < 0.5+g));
			if ( rvecUC[i] < -0.5)
				rvecUC[i] = -0.5;
			if ( rvecUC[i] > 0.5 )
				rvecUC[i] = 0.5;
			//Convert to a coordinate index in the range [0,1[ and fetch the grid index
			rvecUC[i] -= std::floor(rvecUC[i]);
			assert( (rvecUC[i] >= 0) && (rvecUC[i] < 1));
			rvecUC[i] *= primitiveCellGrid.get_grid_dim()[i];
			xyz[i] = std::floor(rvecUC[i]+0.5);
			if ( std::abs(rvecUC[i]-xyz[i]) > 0.01 )
				throw std::logic_error("Could not establish connection between grids of the primitive- and the supercell");
		}
		int cnsq = primitiveCellGrid.get_xyz_to_reducible(xyz);
		assert( (cnsq >= 0) and ( cnsq < primitiveCellGrid.get_num_points() ) );
		rSuperCellToPrimitve[irSC] = std::move(std::make_pair(cnsq, std::move(R) ) );
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

} /* namespace PhononStructure */
} /* namespace elephon */
