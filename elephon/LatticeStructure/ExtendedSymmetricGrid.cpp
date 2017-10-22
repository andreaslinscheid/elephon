/*	This file ExtendedSymmetricGrid.cpp is part of elephon.
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
 *  Created on: Oct 21, 2017
 *      Author: A. Linscheid
 */

#include "LatticeStructure/ExtendedSymmetricGrid.h"
#include <cmath>
#include <set>
#include <map>

namespace elephon
{
namespace LatticeStructure
{

void
ExtendedSymmetricGrid::initialize(
		std::vector<int> dim,
		double gridPrec,
		std::vector<double> shift,
		Symmetry symmetry,
		LatticeModule lattice )
{
	assert( dim.size() == 3 );
	assert( (dim[0] > 1) and ( dim[1] > 1) and (dim[2] > 1)  );
	std::vector<double> extendedUC{ 1.0 / static_cast<double>(dim[0]),
									1.0 / static_cast<double>(dim[1]),
									1.0 / static_cast<double>(dim[2])};

	symmetry_ = std::move( symmetry );
	RegularBareGrid::initialize(
			std::move(dim),
			symmetry_.is_reci(),
			gridPrec,
			std::move(shift),
			std::move(lattice));

	std::set< GridPoint > reducibleSet;
	for ( int i = 0 ; i < this->get_num_points(); ++i)
		reducibleSet.insert( GridPoint( this->get_vector_direct(i), extendedUC, this->get_grid_prec() ));

	// We cannot use SymmetryReduction for finding the irreducible set of grid vectors
	// because this method establishes a set of unique irreducible objects with their stars such that
	// it is possible to retrieve for a give irreducible vector all star vectors which forms the complete
	// grid. Here, the initial set is overcomplete, since vectors (1,0,0) and (0,0,0) are the same by
	// lattice periodicity. All we can achieve is a map reducible to irreducible
	// Since the algorithm thus does not establish a bijective map, it is acceptable that the grid is
	// not closed under its symmetry operations and that stars of k points overlap.
	std::map<GridPoint, std::vector<std::pair<int,int>> > irreducibleSet;
	int ikRred = 0;
	for ( auto redGridVec : reducibleSet)
	{
		bool isSymmetryEquivalentToIrreducible = false;
		int isymOp = symmetry_.get_identity_index();
		for ( int isym = 0 ; isym < symmetry_.get_num_symmetries() ; ++isym )
		{
			auto irredGridVec = redGridVec;
			irredGridVec.transform(symmetry_.get_sym_op(isym));
			auto it = irreducibleSet.find(irredGridVec);
			if ( it != irreducibleSet.end() )
			{
				isSymmetryEquivalentToIrreducible = true;
				isymOp = isym;
				it->second.push_back( std::make_pair(ikRred,isymOp) );
				break;
			}
		}

		if ( ! isSymmetryEquivalentToIrreducible )
		{
			auto ret = irreducibleSet.insert( std::make_pair(redGridVec, std::vector<std::pair<int,int>>() ) );
			ret.first->second.push_back( std::make_pair(ikRred,isymOp) );
		}
		ikRred++;
	}

	nIrredPts_ = irreducibleSet.size();
	irreducibleGridVectorIndices_.reserve(nIrredPts_);
	redToIrredMap_.assign(ikRred, -1);
	symRedToIrredMap_.assign(ikRred, -1);
	for ( auto const & irredVec : irreducibleSet )
	{
		int ikIrreducible = irreducibleGridVectorIndices_.size();
		int ikReducible = irredVec.second.front().first;
		irreducibleGridVectorIndices_.push_back( ikReducible );
		for ( auto ikRedSymPair : irredVec.second )
		{
			redToIrredMap_[ikRedSymPair.first] = ikIrreducible;
			symRedToIrredMap_[ikRedSymPair.first] = ikRedSymPair.second;
		}
	}
	assert( *std::min_element(redToIrredMap_.begin(), redToIrredMap_.end()) == 0 );
	assert( *std::min_element(symRedToIrredMap_.begin(), symRedToIrredMap_.end()) == 0 );
}

std::vector<double>
ExtendedSymmetricGrid::get_vector_direct(int reducibleIndex) const
{
	auto xyz = RegularSymmetricGrid::get_reducible_to_xyz(reducibleIndex);
	auto d = RegularSymmetricGrid::get_grid_dim();
	auto s = RegularSymmetricGrid::get_grid_shift();
	return std::vector<double>{
		(xyz[0] + s[0])/double(d[0]-1) ,
		(xyz[1] + s[1])/double(d[1]-1) ,
		(xyz[2] + s[2])/double(d[2]-1) };
}

namespace detail
{
GridPoint::GridPoint()
{

}

GridPoint::GridPoint(std::vector<double> coords, std::vector<double> extended, double gridPrec)
{
	this->initialize(std::move(coords), std::move(extended), gridPrec);
}

void
GridPoint::initialize(std::vector<double> coords, std::vector<double> extended, double gridPrec)
{
	coords_ = std::move(coords);
	extended_ = std::move(extended);
	gridPrec_ = gridPrec;
	assert(coords_.size() == 3);
	assert(extended_.size() == 3);
	assert( (extended_[0] <= 0.5) and (extended_[1] <= 0.5) and (extended_[2] <= 0.5) );
	this->map_back();
}

void
GridPoint::transform( Symmetry::Sop const & sop )
{
	sop.apply( coords_ , false );
	this->map_back();
}

bool operator< (GridPoint const & d1, GridPoint const & d2)
{
	//Compare position
	assert( d1.gridPrec_ == d2.gridPrec_ );
	for ( int i = 3 ; i--; )
		if ( std::abs( d1.coords_[i] - d2.coords_[i]) > d1.gridPrec_ )
			return d1.coords_[i] < d2.coords_[i];
	return false;
}

void
GridPoint::map_back()
{
	for (int i = 0 ; i < 3 ; ++i )
	{
		double effectiveLength = (1.0 + extended_[i]);
		while ( coords_[i] < 0 )
		{
			coords_[i] += 1.0;
		};
		while ( coords_[i] >= effectiveLength )
		{
			coords_[i] -= 1.0;
		};
	}

}

} /* namespace detail */

} /* namespace LatticeStructure */
} /* namespace elephon */
