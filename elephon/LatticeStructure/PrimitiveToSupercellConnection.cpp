/*	This file PrimitiveToSupercellConnection.cpp is part of elephon.
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
 *  Created on: Feb 2, 2018
 *      Author: A. Linscheid
 */

#include "PrimitiveToSupercellConnection.h"
#include "LatticeStructure/UnitCell.h"
#include "Algorithms/LinearAlgebraInterface.h"
#include "Algorithms/helperfunctions.hpp"
#include "Auxillary/AlignedVector.h"
#include "LatticeStructure/RegularBareGrid.h"
#include <set>

namespace elephon
{
namespace LatticeStructure
{

template<typename T>
struct cmp_vec
{
	cmp_vec(double latPrec) : latPrec_(latPrec) {};

	bool operator() (std::array<T,3> const & a, std::array<T,3> const & b) const
	{
		for ( int i = 3; i-- ; )
			if ( std::abs(a[i]-b[i]) >= latPrec_ )
				return a[i] < b[i];
		return false;
	};

	double latPrec_;
};

void
PrimitiveToSupercellConnection::initialize(
		std::shared_ptr<const UnitCell> primitiveCell,
		std::shared_ptr<const UnitCell> superCell)
{
	assert(primitiveCell->get_atoms_list().size() > 0);
	primitiveCell_ = primitiveCell;
	superCell_ = superCell;
	assert(superCell->get_atoms_list().size() > 0);
	std::fill_n(supercellToPrimitiveCoordsMatrix_.data(), supercellToPrimitiveCoordsMatrix_.size(), 0);
	std::fill_n(primitiveToSupercellCoordsMatrix_.data(), primitiveToSupercellCoordsMatrix_.size(), 0.0);

	auto norm = [] (std::vector<double> const & a) {return std::sqrt(std::pow(a[0],2)+std::pow(a[1],2)+std::pow(a[2],2));};
	int scaleX = std::floor( norm(superCell->get_lattice().get_lattice_vector(0))*superCell->get_lattice().get_alat()
					/(norm(primitiveCell->get_lattice().get_lattice_vector(0))*primitiveCell->get_lattice().get_alat()) + 0.5);
	int scaleY = std::floor( norm(superCell->get_lattice().get_lattice_vector(1))*superCell->get_lattice().get_alat()
					/(norm(primitiveCell->get_lattice().get_lattice_vector(1))*primitiveCell->get_lattice().get_alat()) + 0.5);
	int scaleZ = std::floor( norm(superCell->get_lattice().get_lattice_vector(2))*superCell->get_lattice().get_alat()
					/(norm(primitiveCell->get_lattice().get_lattice_vector(2))*primitiveCell->get_lattice().get_alat()) + 0.5);
	supercellToPrimitiveCoordsMatrix_[0][0] = scaleX;
	supercellToPrimitiveCoordsMatrix_[1][1] = scaleY;
	supercellToPrimitiveCoordsMatrix_[2][2] = scaleZ;

	std::copy(supercellToPrimitiveCoordsMatrix_.data(),
			supercellToPrimitiveCoordsMatrix_.data() + supercellToPrimitiveCoordsMatrix_.size(),
			primitiveToSupercellCoordsMatrix_.data() );

	Algorithms::LinearAlgebraInterface linalg;
	linalg.inverse(primitiveToSupercellCoordsMatrix_, primitiveToSupercellCoordsMatrix_);

	for (int ia = 0 ; ia < primitiveCell_->get_atoms_list().size(); ++ia)
		lookupPrimitive_.insert(std::make_pair(primitiveCell_->get_atoms_list()[ia], ia));

	for (int ia = 0 ; ia < superCell->get_atoms_list().size(); ++ia)
		lookupSupercell_.insert(std::make_pair(superCell->get_atoms_list()[ia], ia));

	// to make the connection we have to work in cartesian coordinates. To make search efficient,
	// we define a comparison function similar to the one of LatticeStructure::Atom
	const double latPrec = primitiveCell_->get_atoms_list().front().get_position_precision();

	indexRVectorMap_.resize(boost::extents[scaleX][scaleY][scaleZ]);
	std::fill_n(indexRVectorMap_.data(), indexRVectorMap_.size(), -1);
	std::set<std::array<int,3>, cmp_vec<int>> rvectors(latPrec);
	auto hint = rvectors.end();
	for (int iRx = 0; iRx < scaleX; ++iRx)
		for (int iRy = 0; iRy < scaleY; ++iRy)
			for (int iRz = 0; iRz < scaleZ; ++iRz)
			{
				hint = rvectors.insert(hint, std::array<int,3>({iRx,iRy,iRz}));
			}
	Rvectors_.resize(boost::extents[rvectors.size()][3]);
	int iR = 0;
	for (auto & R : rvectors)
	{
		assert(R.size() == 3);
		for (int ix = 0 ; ix < 3; ++ix)
			Rvectors_[iR][ix] = R[ix];
		indexRVectorMap_[R[0]][R[1]][R[2]] = iR;
		++iR;
	}

	superCellToPrimitveAtomIndex_.resize(superCell_->get_atoms_list().size());
	superCellToPrimitveRVector_.resize(boost::extents[superCell_->get_atoms_list().size()][3]);
	for (int iaSC = 0 ; iaSC < superCell_->get_atoms_list().size(); ++iaSC)
	{
		std::array<int,3> R;
		this->compute_get_equiv_atom_primitive(iaSC, R, superCellToPrimitveAtomIndex_[iaSC]);
		for (int xi = 0 ; xi < 3; ++xi)
			superCellToPrimitveRVector_[iaSC][xi] = R[xi];
	}
}

void
PrimitiveToSupercellConnection::discover_primitive_cell(
		std::shared_ptr<const UnitCell> superCell,
		UnitCell & primitiveCell) const
{
	// a primitive cell contained in a supercell is found by computing all
	// vectors from a given atom to all atoms of the same type within the cell.
	// If any of those vectors maps the cell into itself up to periodic boundary conditions,
	// we have found a lattice vector for the reduced cell.
}

void
PrimitiveToSupercellConnection::get_supercell_vectors(
		Auxillary::Multi_array<int,2> & Rvectors) const
{
	assert(Rvectors_.size() > 0);
	Rvectors.resize(boost::extents[Rvectors_.shape()[0]][3]);
	std::copy(Rvectors_.data(), Rvectors_.data()+Rvectors_.size(), Rvectors.data());
}

std::array<int,3>
PrimitiveToSupercellConnection::get_supercell_vector(int iR) const
{
	assert((iR>=0)&&(iR<Rvectors_.size()));
	return std::array<int,3>({Rvectors_[iR][0], Rvectors_[iR][1], Rvectors_[iR][2]});
}

int
PrimitiveToSupercellConnection::get_equiv_atom_primitive(int iASC) const
{
	assert((iASC>=0)&&(iASC<superCellToPrimitveAtomIndex_.size()));
	return superCellToPrimitveAtomIndex_[iASC];
}

int
PrimitiveToSupercellConnection::get_equiv_atom_primitive_lattice_vector_index(int iASC) const
{
	assert((iASC>=0)&&(iASC<superCellToPrimitveRVector_.shape()[0]));
	const int Rx = superCellToPrimitveRVector_[iASC][0];
	const int Ry = superCellToPrimitveRVector_[iASC][1];
	const int Rz = superCellToPrimitveRVector_[iASC][2];
	const int iR = indexRVectorMap_[Rx][Ry][Rz];
	return iR;
}

void
PrimitiveToSupercellConnection::compute_get_equiv_atom_primitive(int iASC, std::array<int,3> & R, int & iAPC) const
{
	assert((superCell_->get_atoms_list().size()>iASC) && (iASC>=0));
	auto atom = superCell_->get_atoms_list()[iASC];
	auto pos = atom.get_position();
	this->supercell_to_primitive_coordinates(pos);
	// map back to 1. primitive cell
	for (int i = 0 ; i < 3; ++i)
	{
		R[i] = std::floor(pos[i]+0.5);
		pos[i] -= R[i] ;
	}
	atom.set_position(std::move(pos));
	iAPC = this->find_atom(atom,true);
}

int
PrimitiveToSupercellConnection::find_atom(Atom const & a, bool primitiveCell) const
{
	if ( primitiveCell )
	{
		auto it = lookupPrimitive_.find(a);
		if ( it == lookupPrimitive_.end())
			return -1;
		return it->second;
	}
	else
	{
		auto it = lookupSupercell_.find(a);
		if ( it == lookupSupercell_.end())
			return -1;
		return it->second;
	}
}

void
PrimitiveToSupercellConnection::supercell_to_primitive_coordinates(std::vector<double> & vec) const
{
	this->supercell_to_primitive_coordinates_no_shift(vec);
}

void
PrimitiveToSupercellConnection::get_embedded_supercell_vectors(
		Auxillary::Multi_array<int,2> & Rvectors,
		std::array<int,3> & MaxEachDirSpace,
		Auxillary::Multi_array<int,3> & MapRToIndex) const
{
	// for the currently implemented functionality, this is fairly easy.
	this->get_supercell_vectors(Rvectors);

	std::array<std::pair<int,int>, 3> dims;
	for (int iR = 0; iR < Rvectors.shape()[0] ; ++iR)
		for (int ix = 0; ix < 3; ++ix)
		{
			dims[ix].first = std::min(dims[ix].first, Rvectors[iR][ix]);
			dims[ix].second = std::max(dims[ix].second, Rvectors[iR][ix]);
		}

	// shift to the center
	for (int ix = 0; ix < 3; ++ix)
	{
		MaxEachDirSpace[ix] = (dims[ix].second - dims[ix].first)/2 + 1;
		for (int iR = 0; iR < Rvectors.shape()[0] ; ++iR)
		{
			Rvectors[iR][ix] -= (dims[ix].second + dims[ix].first)/2;
			assert((-MaxEachDirSpace[ix] < Rvectors[iR][ix]) && (Rvectors[iR][ix] < MaxEachDirSpace[ix]));
		}
	}

	// if we have unbalanced R vectors, add the inverse to the set and finally copy
	// the balanced lattice vectors to Rvectors
	const double latPrec = primitiveCell_->get_atoms_list().front().get_position_precision();
	std::set<std::array<int,3>, cmp_vec<int> > balancedRset(latPrec);
	for (int iRUB = 0 ;iRUB < Rvectors.shape()[0]; ++iRUB)
	{
		balancedRset.insert(std::array<int,3>({ Rvectors[iRUB][0],  Rvectors[iRUB][1],  Rvectors[iRUB][2]}));
		balancedRset.insert(std::array<int,3>({-Rvectors[iRUB][0], -Rvectors[iRUB][1], -Rvectors[iRUB][2]}));
	}
	Rvectors.resize(boost::extents[balancedRset.size()][3]);
	typedef boost::multi_array_types::extent_range range;
	MapRToIndex.resize( boost::extents[range(-MaxEachDirSpace[0],MaxEachDirSpace[0])]
									  [range(-MaxEachDirSpace[1],MaxEachDirSpace[1])]
									  [range(-MaxEachDirSpace[2],MaxEachDirSpace[2])]);
	std::fill_n(MapRToIndex.data(), MapRToIndex.size(), -1);
	auto it = balancedRset.begin();
	for (int iRB = 0 ;iRB < balancedRset.size(); ++iRB)
	{
		for (int ix = 0; ix < 3; ++ix)
			Rvectors[iRB][ix] = (*it)[ix];
		MapRToIndex[(*it)[0]][(*it)[1]][(*it)[2]] = iRB;
		++it;
	}
}

void
PrimitiveToSupercellConnection::primitive_to_supercell_coordinates(std::vector<double> & vec) const
{
	this->primitive_to_supercell_coordinates_no_shift(vec);
}

int
PrimitiveToSupercellConnection::primitive_to_supercell_atom_index(int primitiveAtomIndex) const
{
	assert((primitiveAtomIndex>=0)&&(primitiveAtomIndex<primitiveCell_->get_atoms_list().size()));
	auto a = primitiveCell_->get_atoms_list()[primitiveAtomIndex];
	auto pos = a.get_position();
	this->primitive_to_supercell_coordinates(pos);
	a.set_position(pos);
	int indexSC = this->find_atom(a, false);
	assert(indexSC>=0);
	return indexSC;
}

int
PrimitiveToSupercellConnection::get_supercell_volume_factor() const
{
	int fact = Algorithms::helperfunctions::determinant_3by3_matrix(supercellToPrimitiveCoordsMatrix_);
	return fact;
}

void
PrimitiveToSupercellConnection::build_supercell_to_primitive(
		LatticeStructure::RegularBareGrid const & primitiveCellGrid,
		LatticeStructure::RegularBareGrid const & supercellGrid,
		std::vector<LatticeStructure::Atom> const & atomsUC,
		Auxillary::Multi_array<int,2> & indexMap,
		Auxillary::Multi_array<int,3> & latticeVectorMap,
		std::array<std::pair<int,int>,3> & Rdims,
		Auxillary::Multi_array<int,2> & listAllRVectors,
		Auxillary::Multi_array<int,3> & RVectorIndexMap ) const
{
	const int nRSC = supercellGrid.get_num_points();
	const double gp = primitiveCellGrid.get_grid_prec();

	auto rVectorsSupercell = supercellGrid.get_all_vectors_grid();

	const int nA = static_cast<int>(atomsUC.size());
	indexMap.resize(boost::extents[nA][nRSC]);
	latticeVectorMap.resize(boost::extents[nA][nRSC][3]);
	std::vector<int> xyz(3);
	decltype(rVectorsSupercell) rVectorsShiftedSupercell;
	for ( int ia = 0 ; ia < atomsUC.size() ; ++ia )
	{
		// transform to a coordinate system, where the current atom is in the center and apply periodic
		// boundary conditions.
		rVectorsShiftedSupercell.assign(rVectorsSupercell.begin(), rVectorsSupercell.end());
		auto aPos = atomsUC[ia].get_position();
		this->primitive_to_supercell_coordinates_no_shift(aPos);
		for ( int irSC = 0 ; irSC < nRSC ; ++irSC )
			for ( int i = 0 ; i < 3 ; ++i )
			{
				rVectorsShiftedSupercell[irSC*3+i] -= aPos[i];
				rVectorsShiftedSupercell[irSC*3+i] -= std::floor(rVectorsShiftedSupercell[irSC*3+i]+0.5);
			}

		// go to coordinates of the primitive cell
		this->supercell_to_primitive_coordinates_no_shift(rVectorsShiftedSupercell);

		// in this coordinates, pick the first unit cell, which is now of cause centered around
		// the particular atom atomsUC[ia]
		for ( int irSC = 0 ; irSC < nRSC ; ++irSC )
		{
			for ( int i = 0 ; i < 3 ; ++i )
			{
				latticeVectorMap[ia][irSC][i] = std::floor(rVectorsShiftedSupercell[irSC*3+i]+0.5);
				double rPC = rVectorsShiftedSupercell[irSC*3+i] - latticeVectorMap[ia][irSC][i];
				// keep track of the minima and maxima of the lattice vectors in each direction
				Rdims[i].first = std::min(Rdims[i].first, latticeVectorMap[ia][irSC][i]);
				Rdims[i].second = std::max(Rdims[i].second, latticeVectorMap[ia][irSC][i]);

				// numerical inaccuracy may occur. For margins within the grid precision,
				// we take the liberty to reset so that the resulting vectors is strictly in
				// the first Brillouin zone (but the right border included!).
				assert( (rPC >= -0.5-gp)
						&& (rPC < 0.5+gp));
				if ( rPC < -0.5)
					rPC = -0.5;
				if ( rPC > 0.5 )
					rPC = 0.5;

				//Convert to a coordinate index in the range [0,1] and fetch the grid index
				rPC -= std::floor(rPC);
				rPC *= primitiveCellGrid.get_grid_dim()[i];
				xyz[i] = std::floor(rPC+0.5);
			}
			// NOTE: since the atom position shift means we are out of the grid, we have to allow for
			// the case where any xyz[i] is equal to primitiveCellGrid.get_grid_dim()[i]. Thus we have to use the
			// periodic version get_xyz_to_reducible.
			int cnsq = primitiveCellGrid.get_xyz_to_reducible_periodic(xyz);
			assert( (cnsq >= 0) and ( cnsq < primitiveCellGrid.get_num_points() ) );
			indexMap[ia][irSC] = cnsq;
		}
	}
}

Auxillary::Multi_array<int,2> const &
PrimitiveToSupercellConnection::get_supercell_matrix() const
{
	return superCellToPrimitveRVector_;
}

void
PrimitiveToSupercellConnection::supercell_to_primitive_coordinates_no_shift(std::vector<double> & vec) const
{
	assert(vec.size() % 3 == 0);
	for (int iv = 0; iv < vec.size()/3 ; ++iv)
	{
		double buff[] = {vec[iv*3+0], vec[iv*3+1], vec[iv*3+2]};
		for (int i = 0 ; i < 3; ++i)
			vec[iv*3+i] =supercellToPrimitiveCoordsMatrix_[i][0]*buff[0]+
						 supercellToPrimitiveCoordsMatrix_[i][1]*buff[1]+
						 supercellToPrimitiveCoordsMatrix_[i][2]*buff[2];
	}
}

void
PrimitiveToSupercellConnection::primitive_to_supercell_coordinates_no_shift(std::vector<double> & vec) const
{
	assert(vec.size() % 3 == 0);
	for (int iv = 0; iv < vec.size()/3 ; ++iv)
	{
		double buff[] = {vec[iv*3+0], vec[iv*3+1], vec[iv*3+2]};
		for (int i = 0 ; i < 3; ++i)
			vec[iv*3+i] =primitiveToSupercellCoordsMatrix_[i][0]*buff[0]+
						 primitiveToSupercellCoordsMatrix_[i][1]*buff[1]+
						 primitiveToSupercellCoordsMatrix_[i][2]*buff[2];
	}
}

} /* namespace LatticeStructure */
} /* namespace elephon */
