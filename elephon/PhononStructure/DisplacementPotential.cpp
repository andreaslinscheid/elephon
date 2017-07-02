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

namespace elephon
{
namespace PhononStructure
{

void
DisplacementPotential::build(  LatticeStructure::UnitCell const & unitCell,
		LatticeStructure::UnitCell const & superCell,
		std::vector<LatticeStructure::AtomDisplacement> const & irredDispl,
		LatticeStructure::RegularGrid const & unitcellGrid,
		LatticeStructure::RegularGrid const & supercellGrid,
		std::vector<double> const & potentialUC,
		std::vector< std::vector<double> > const & potentialDispl )
{
	unitCell.compute_supercell_dim(superCell,superCellDim_);
	numModes_ = unitCell.get_atoms_list().size()*3;
	int nRSC = supercellGrid.get_np_red();
	int nR = nRSC/this->get_num_R();
	assert( unitcellGrid.get_np_red() == nR );
	assert( potentialDispl.size() == irredDispl.size() );

	//Create a table which maps a point in real space in the supercell
	// into a point in the primitive cell plus a lattice vector
	std::vector< std::pair<int,std::vector<int> > > rSuperCellToPrimitve( nRSC );
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
			rvecUC[i] = (rvec[i]-R[i])*unitcellGrid.get_grid_dim()[i];
			xyz[i] = std::floor(rvecUC[i]+0.5);
			if ( std::abs(rvecUC[i]-xyz[i]) > 0.01 )
				throw std::logic_error("Could not establish connection between grids of the primitive- and the supercell");
		}
		int cnsq = unitcellGrid.get_xyz_to_reducible(xyz);
		assert( cnsq < unitcellGrid.get_np_red() );
		rSuperCellToPrimitve[irSC] = std::move(std::make_pair(cnsq, std::move(R) ) );
	}

	//From the potentials, construct the difference
	std::vector< std::vector<double> > potentialVariation = potentialDispl;
	for (int ird = 0 ; ird < potentialDispl.size(); ++ird)
		for (int ir = 0 ; ir < nRSC; ++ir)
			potentialVariation[ird][ir] -= potentialUC[rSuperCellToPrimitve[ir].first];

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
					symmetryShiftedGrid[ir*3+j] *= supercellGrid.get_grid_dim()[j];
					xyz[j] = std::floor(symmetryShiftedGrid[ir*3+j]+0.5);
					if ( std::abs(symmetryShiftedGrid[ir*3+j]-xyz[j]) > 0.01 )
						throw std::logic_error("The symmetry operation that transports the atom "
								"does not map the grid to itself which can't be.");
				}
				int cnsq = supercellGrid.get_xyz_to_reducible(xyz);
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
}

void
DisplacementPotential::compute_rot_map(
		std::vector<double> const & shift,
		LatticeStructure::RegularGrid const & supercellGrid,
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
				shiftedGrid[ir*3+j] *= supercellGrid.get_grid_dim()[j];
				xyz[j] = std::floor(shiftedGrid[ir*3+j]+0.5);
				if ( std::abs(shiftedGrid[ir*3+j]-xyz[j]) > 0.01 )
					throw std::logic_error("The symmetry operation does not map the grid to itself which can't be.");
			}
			int cnsq = supercellGrid.get_xyz_to_reducible(xyz);
			rotMap[isym][cnsq] = ir;
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

} /* namespace PhononStructure */
} /* namespace elephon */
