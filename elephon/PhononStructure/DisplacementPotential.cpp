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

namespace elephon
{
namespace PhononStructure
{

void
DisplacementPotential::build(  LatticeStructure::UnitCell unitCell,
		LatticeStructure::UnitCell const & superCell,
		std::vector<LatticeStructure::AtomDisplacement> const & irredDispl,
		LatticeStructure::RegularSymmetricGrid unitcellGrid,
		LatticeStructure::RegularSymmetricGrid const & supercellGrid,
		std::vector<double> const & potentialUC,
		std::vector< std::vector<double> > const & potentialDispl )
{
	unitCell.compute_supercell_dim(superCell,superCellDim_);
	numModes_ = unitCell.get_atoms_list().size()*3;
	int nRSC = supercellGrid.get_np_red();
	int nR = nRSC/this->get_num_R();
	nptsRealSpace_ = nR;
	assert( unitcellGrid.get_np_red() == nR );
	assert( potentialDispl.size() == irredDispl.size() );

	std::vector< std::pair<int,std::vector<int> > > rSuperCellToPrimitve;
	this->build_supercell_to_primite(unitcellGrid, supercellGrid, rSuperCellToPrimitve);
	//Create a table which maps a point in real space in the supercell
	// into a point in the primitive cell plus a lattice vector
	for ( int irSC = 0 ; irSC < nRSC ; ++irSC )
	{
		auto rvec = supercellGrid.get_vector_direct(irSC);
		std::vector<int> R(3);
		std::vector<double> rvecUC(3);
		std::vector<int> xyz(3);
		for ( int i = 0 ; i < 3 ; ++i )
		{
			//scale to units of the primitive cell
			rvec[i] *= superCellDim_[i];
			//We add a tiny bit to the vector so that 1.9999999 will not (incorrectly) map
			//to R = 1 and then 0.9999999 will (correctly) map to the lattice site at 0.
			R[i] = std::floor(rvec[i]+0.5+2*unitCell.get_symmetry().get_symmetry_prec());
			rvecUC[i] = rvec[i]-R[i];
			assert( (rvecUC[i] >= -0.5) && (rvecUC[i] < 0.5));
			//Convert to a coordinate index in the range [0,1[
			rvecUC[i] -= std::floor(rvecUC[i]);
			assert( (rvecUC[i] >= 0) && (rvecUC[i] < 1));
			rvecUC[i] *= unitcellGrid.get_grid_dim()[i];
			xyz[i] = std::floor(rvecUC[i]+0.5);
			if ( std::abs(rvecUC[i]-xyz[i]) > 0.01 )
				throw std::logic_error("Could not establish connection between grids of the primitive- and the supercell");
		}
		int cnsq = unitcellGrid.get_xyz_to_reducible(xyz);
		assert( (cnsq >= 0) and ( cnsq < unitcellGrid.get_np_red() ) );
		rSuperCellToPrimitve[irSC] = std::move(std::make_pair(cnsq, std::move(R) ) );
	}

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
			unitCell.get_symmetry(),
			unitCell.get_atoms_list(),  irredAtoms,
			redToIrredAtoms, symRedToIrredAtoms,
			irredToRedAtoms, symIrredToRedAtoms);

	//define a consistent atom numbering in the primitive cell.
	//We choose the unitCell.get_atoms_list() order
	std::map< LatticeStructure::Atom, int > primitiveCellLookup;
	for ( int ia = 0 ; ia < unitCell.get_atoms_list().size(); ++ia )
		primitiveCellLookup.insert( std::move( std::make_pair( unitCell.get_atoms_list()[ia], ia ) ) );

	int naSC = superCell.get_atoms_list().size();

	Algorithms::LinearAlgebraInterface linAlg;

	//For each such atom, compute the local displacement potential in x,y and z
	std::vector<double> irreducibleDisplPot( irredAtoms.size()*3*nRSC );

	//This counter advances the point in the set of input potentialVariations.
	//At each irreducible atomic site, we add the irreducible displacements.
	int iIrredDisplCount = 0;
	for ( int ia = 0 ; ia < irredAtoms.size(); ++ia )
	{
		//regenerate displacements to make the connection with the input forces
		LatticeStructure::Symmetry const & sSym = unitCell.get_site_symmetry(ia);

		std::vector<LatticeStructure::AtomDisplacement> irredNotUsed,reducible;
		bool symmetricDispl = not irredDispl[iIrredDisplCount].is_plus_minus_displ_equivalent();
		std::vector<int> redToIrredDispl, symRedToIrredDispl;
		std::vector< std::vector<int> > irredToRedDispl, symIrredToRedDispl;
		unitCell.get_site_displacements( irredAtoms[ia],
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
			auto d = irredDispl[iIrredDisplCount+idir].get_direction();
			for ( auto &xi : d )
				xi *= irredDispl[iIrredDisplCount].get_magnitude();
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
		for ( int i = 0 ; i < 3*nRSC; ++i)
			linDVscf[i] = irreducibleDisplPot[irA*3*nRSC+i];

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
			superCell.get_symmetry().apply(isym,symmetryShiftedGrid.begin(),symmetryShiftedGrid.end(),true);
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
				for (int xi = 0 ; xi < 3 ; ++xi)
					linDVscfTransform[cnsq*3+xi]=linDVscf[ir*3+xi];
			}
			superCell.get_symmetry().rotate_cartesian(
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
	unitCell_ = std::move(unitCell);
	unitCellGrid_ = std::move(unitcellGrid);
}

void
DisplacementPotential::compute_rot_map(
		std::vector<double> const & shift,
		LatticeStructure::RegularSymmetricGrid const & supercellGrid,
		LatticeStructure::Symmetry const & siteSymmetry,
		std::vector< std::vector<int> > & rotMap) const
{
	int nR = supercellGrid.get_np_red();
	int iS = siteSymmetry.get_num_symmetries();
	rotMap = std::vector< std::vector<int> >(iS, std::vector<int>(nR) );
	std::vector<double> shiftedGrid(nR*3);
	for ( int i = 0 ; i < nR ; ++i)
	{
		auto r = supercellGrid.get_vector_direct(i);
		for ( int xi = 0 ; xi < 3; ++xi)
			shiftedGrid[i*3+xi] = r[xi] - shift[xi];
	}

	for ( int isym = 0 ; isym < iS; ++isym)
	{
		//rotate all grid points
		siteSymmetry.rotate<double>(isym,shiftedGrid.begin(),shiftedGrid.end(),true);
		//since the grid is regular, we must obtain the (x,y,z) indices by scaling with the grid dimension
		//in each direction
		std::vector<int> xyz(3);
		for ( int ir = 0 ; ir < nR; ++ir)
		{
			for ( int j = 0 ; j < 3; ++j)
			{
				shiftedGrid[ir*3+j] -= std::floor(shiftedGrid[ir*3+j]);
				shiftedGrid[ir*3+j] *= supercellGrid.get_grid_dim()[j];
				xyz[j] = std::floor(shiftedGrid[ir*3+j]+0.5);
				if ( std::abs(shiftedGrid[ir*3+j]-xyz[j]) > 0.01 )
					throw std::logic_error("The symmetry operation does not map the grid to itself which can't be.");

			}
			int cnsq = supercellGrid.get_xyz_to_reducible(xyz);
			assert( (cnsq >= 0) and (cnsq < nR));
			rotMap[isym][cnsq] = ir;
		}
	}
}

void
DisplacementPotential::compute_dvscf_q(
		std::vector<double> const & qVect,
		std::vector<std::complex<double>> const & dynamicalMatrices,
		std::vector<double> const & masses,
		std::vector<std::complex<float>> & dvscf) const
{
	assert( qVect.size() % 3 == 0 );
	int nq = qVect.size()/3;

	Algorithms::LinearAlgebraInterface linAlg;

	assert( (dynamicalMatrices.size()/nq) / (numModes_*numModes_) == 1 );
	assert( (masses.size()*3) / numModes_ == 1 );

	dvscf.resize(nq*numModes_*nptsRealSpace_);
	std::vector<std::complex<float>> ftDisplPot(nq*numModes_*nptsRealSpace_ , std::complex<float>(0) );
	std::vector<std::complex<float> > dynmatMass(numModes_*numModes_);
	for ( int iq = 0 ; iq < nq; ++iq)
	{
		for ( int iRz = 0 ; iRz < superCellDim_[2]; ++iRz )
			for ( int iRy = 0 ; iRy < superCellDim_[1]; ++iRy )
				for ( int iRx = 0 ; iRx < superCellDim_[0]; ++iRx )
				{
					int R[] = {	iRx <= superCellDim_[0]/2 ? iRx : iRx - superCellDim_[0],
									iRy <= superCellDim_[1]/2 ? iRy : iRy - superCellDim_[1],
									iRz <= superCellDim_[2]/2 ? iRz : iRz - superCellDim_[2] };
					std::complex<float> phase = std::complex<float>(std::exp( std::complex<double>(0,
							2.0*M_PI*(qVect[iq*3+0]*R[0]+qVect[iq*3+1]*R[1]+qVect[iq*3+2]*R[2]) )));

					for ( int mu = 0 ; mu < numModes_; ++mu)
						for ( int ir = 0 ; ir < nptsRealSpace_; ++ir)
							ftDisplPot[(iq*numModes_+mu)*nptsRealSpace_+ir]
									   += phase*data_[this->mem_layout(ir,mu,this->RVectorLayout(iRx,iRy,iRz))];
				}

		for ( int mu1 = 0 ; mu1 < numModes_; ++mu1 )
			for ( int mu2 = 0 ; mu2 < numModes_; ++mu2 )
				dynmatMass[mu1*numModes_+mu2] = dynamicalMatrices[(iq*numModes_+mu1)*numModes_+mu2] / masses[mu2/3];

		auto r_ptr = &dvscf[iq*numModes_*nptsRealSpace_];
		auto dmat_ptr = &dynmatMass[0];
		auto ftdp_ptr = &ftDisplPot[iq*numModes_*nptsRealSpace_];
		linAlg.call_gemm( 'n', 'n',
				numModes_, nptsRealSpace_ , numModes_,
				std::complex<float>(1.0f), dmat_ptr, numModes_,
				ftdp_ptr, nptsRealSpace_,
				std::complex<float>(0.0f), r_ptr, nptsRealSpace_);
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
	return ir+nptsRealSpace_*(mu+numModes_+iR);
}

int
DisplacementPotential::RVectorLayout(int iRz, int iRy, int iRx ) const
{
	assert( (iRx >= 0) && (iRx < superCellDim_[0]) );
	assert( (iRy >= 0) && (iRy < superCellDim_[1]) );
	assert( (iRz >= 0) && (iRz < superCellDim_[2]) );
	return (iRz*superCellDim_[1]+iRy)*superCellDim_[0]+iRx;
}

LatticeStructure::RegularSymmetricGrid const &
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
	auto sc = unitCell_.build_supercell( superCellDim_[0], superCellDim_[1], superCellDim_[2] );
	elephon::LatticeStructure::RegularSymmetricGrid scGrid;
	auto dim = unitCellGrid_.get_grid_dim();
	scGrid.initialize(
			std::vector<int>({dim[0]*superCellDim_[0],dim[1]*superCellDim_[1],dim[2]*superCellDim_[2]}),
			unitCellGrid_.get_grid_prec(),
			unitCellGrid_.get_grid_shift(),
			sc.get_symmetry(),
			sc.get_lattice() );

	std::vector< std::pair<int,std::vector<int> > > rSuperCellToPrimitve;
	this->build_supercell_to_primite( unitCellGrid_, scGrid, rSuperCellToPrimitve);

	std::vector<double> dataSC( scGrid.get_np_red() );
	assert( rSuperCellToPrimitve.size() == scGrid.get_np_red() );
	for ( int irSC = 0 ; irSC < scGrid.get_np_red(); ++irSC )
	{
		auto ir = rSuperCellToPrimitve[irSC].first;
		auto & R = rSuperCellToPrimitve[irSC].second;
		int iRx = R[0] < 0 ? R[0] + superCellDim_[0] : R[0];
		int iRy = R[1] < 0 ? R[1] + superCellDim_[1] : R[1];
		int iRz = R[2] < 0 ? R[2] + superCellDim_[2] : R[2];
		dataSC[irSC] = data_[this->mem_layout( ir , atomIndex*3+xi, this->RVectorLayout(iRz,iRy,iRx))];
	}
	writer.write_file(filename, comment, scGrid.get_grid_dim(), sc, dataSC );
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
	this->compute_dvscf_q(qVect,dynamicalMatrices,masses,qDataPlus);

	int nr = unitCellGrid_.get_np_red();
	assert(qDataPlus.size() == nq*numModes_*nr);

	auto mod_filename = [] (std::string filename, int iq, int mu) {
		if ( filename.length() >= 5 )
			if ( filename.substr(filename.length()-4,filename.length()).compare( ".dat" ) == 0 )
				filename.insert(filename.length()-4,"_"+std::to_string(iq)+"_"+std::to_string(mu));
		return filename;
	};

	IOMethods::WriteVASPRealSpaceData writer;
	for ( int iq = 0; iq < qVect.size()/3 ; ++iq )
	{
		for ( auto mu : modeIndices )
		{
			assert( (mu >= 0) && (mu < numModes_) );
			std::string comment = "Displacement potential dvscf(r)/du(q,mu); q = ("
					+std::to_string(qVect[iq*3+0])+" , "+std::to_string(qVect[iq*3+1])+" , "+std::to_string(qVect[iq*3+2])
					+") module # "+ std::to_string(mu) +"\n";
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
		LatticeStructure::RegularSymmetricGrid const & primitiveCellGrid,
		LatticeStructure::RegularSymmetricGrid const & supercellGrid,
		std::vector< std::pair<int,std::vector<int> > > & rSuperCellToPrimitve) const
{
	int nRSC = supercellGrid.get_np_red();

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
			//scale to units of the primitive cell
			rvec[i] *= superCellDim_[i];
			//We add a tiny bit to the vector so that 1.9999999 will not (incorrectly) map
			//to R = 1 and then 0.9999999 will (correctly) map to the lattice site at 0.
			R[i] = std::floor(rvec[i]+0.5+2*primitiveCellGrid.get_grid_prec());
			rvecUC[i] = rvec[i]-R[i];
			assert( (rvecUC[i] >= -0.5) && (rvecUC[i] < 0.5));
			//Convert to a coordinate index in the range [0,1[
			rvecUC[i] -= std::floor(rvecUC[i]);
			assert( (rvecUC[i] >= 0) && (rvecUC[i] < 1));
			rvecUC[i] *= primitiveCellGrid.get_grid_dim()[i];
			xyz[i] = std::floor(rvecUC[i]+0.5);
			if ( std::abs(rvecUC[i]-xyz[i]) > 0.01 )
				throw std::logic_error("Could not establish connection between grids of the primitive- and the supercell");
		}
		int cnsq = primitiveCellGrid.get_xyz_to_reducible(xyz);
		assert( (cnsq >= 0) and ( cnsq < primitiveCellGrid.get_np_red() ) );
		rSuperCellToPrimitve[irSC] = std::move(std::make_pair(cnsq, std::move(R) ) );
	}
}

} /* namespace PhononStructure */
} /* namespace elephon */
