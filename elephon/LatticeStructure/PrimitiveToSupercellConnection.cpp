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
#include <set>

namespace elephon
{
namespace LatticeStructure
{

template<typename T>
struct cmp_vec
{
	cmp_vec(double latPrec) : latPrec_(latPrec) {};

	bool operator() (std::vector<T> const & a, std::vector<T> const & b) const
	{
		assert((a.size()==3)&&(b.size()==3));
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
	primitiveAtomList_ = primitiveCell->get_atoms_list();
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

	for (int ia = 0 ; ia < primitiveAtomList_.size(); ++ia)
		lookupPrimitive_.insert(std::make_pair(primitiveAtomList_[ia], ia));

	for (int ia = 0 ; ia < superCell->get_atoms_list().size(); ++ia)
		lookupSupercell_.insert(std::make_pair(superCell->get_atoms_list()[ia], ia));

	// to make the connection we have to work in cartesian coordinates. To make search efficient,
	// we define a comparison function similar to the one of LatticeStructure::Atom
	const double latPrec = primitiveAtomList_.front().get_position_precision();

	std::set<std::vector<int>, cmp_vec<int>> rvectors(latPrec);
	auto hint = rvectors.end();
	for (int iRx = 0; iRx < scaleX; ++iRx)
		for (int iRy = 0; iRy < scaleY; ++iRy)
			for (int iRz = 0; iRz < scaleZ; ++iRz)
			{
				hint = rvectors.insert(hint, {iRx,iRy,iRz});
			}
	Rvectors_.resize(boost::extents[rvectors.size()][3]);
	int iR = 0;
	for (auto & R : rvectors)
	{
		assert(R.size() == 3);
		for (int ix = 0 ; ix < 3; ++ix)
			Rvectors_[iR][ix] = R[ix];
		++iR;
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
PrimitiveToSupercellConnection::supercell_vectors(
		Auxillary::Multi_array<int,2> & Rvectors) const
{
	assert(Rvectors_.size() > 0);
	Rvectors.resize(boost::extents[Rvectors_.shape()[0]][3]);
	std::copy(Rvectors_.data(), Rvectors_.data()+Rvectors_.size(), Rvectors.data());
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
	assert(vec.size() == 3);
	double buff[] = {vec[0], vec[1], vec[2]};
	for (int i = 0 ; i < 3; ++i)
		vec[i] =supercellToPrimitiveCoordsMatrix_[i][0]*buff[0]+
				supercellToPrimitiveCoordsMatrix_[i][1]*buff[1]+
				supercellToPrimitiveCoordsMatrix_[i][2]*buff[2];
}

void
PrimitiveToSupercellConnection::primitive_to_supercell_coordinates(std::vector<double> & vec) const
{
	assert(vec.size() == 3);
	double buff[] = {vec[0], vec[1], vec[2]};
	for (int i = 0 ; i < 3; ++i)
		vec[i] =primitiveToSupercellCoordsMatrix_[i][0]*buff[0]+
				primitiveToSupercellCoordsMatrix_[i][1]*buff[1]+
				primitiveToSupercellCoordsMatrix_[i][2]*buff[2];
}

int
PrimitiveToSupercellConnection::primitive_to_supercell_atom_index(int primitiveAtomIndex) const
{
	assert((primitiveAtomIndex>=0)&&(primitiveAtomIndex<primitiveAtomList_.size()));
	auto a = primitiveAtomList_[primitiveAtomIndex];
	auto pos = a.get_position();
	this->primitive_to_supercell_coordinates(pos);
	a.set_position(pos);
	int indexSC = this->find_atom(a, false);
	assert(indexSC>=0);
	return indexSC;
}

int
PrimitiveToSupercellConnection::supercell_volume_factor() const
{
	int fact = Algorithms::helperfunctions::determinant_3by3_matrix(supercellToPrimitiveCoordsMatrix_);
	return fact;
}


} /* namespace LatticeStructure */
} /* namespace elephon */
