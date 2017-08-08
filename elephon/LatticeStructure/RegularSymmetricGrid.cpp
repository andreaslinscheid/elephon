/*	This file RegularGrid.cpp is part of elephon.
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
 *  Created on: May 20, 2017
 *      Author: A. Linscheid
 */

#include "SymmetryReduction.h"
#include <assert.h>
#include <LatticeStructure/RegularSymmetricGrid.h>
#include <cmath>
#include <stdexcept>
#include <iostream>

namespace elephon
{
namespace LatticeStructure
{

RegularBareGrid
RegularSymmetricGrid::view_bare_grid() const
{
	return *static_cast<RegularBareGrid const *>(this);
}

void
RegularSymmetricGrid::initialize(
		std::vector<int> dim,
		double gridPrec,
		std::vector<double> shift,
		Symmetry symmetry,
		LatticeModule lattice,
		std::vector<double> const & constraintIrreducible)
{
	symmetry_ = std::move( symmetry );
	RegularBareGrid::initialize(
			std::move(dim),
			symmetry_.is_reci(),
			gridPrec,
			std::move(shift),
			std::move(lattice));

	std::map< RegularSymmetricGrid::GridPoint, int > reducibleSet;
	this->construct_grid_vector_set( reducibleSet );
	std::vector< RegularSymmetricGrid::GridPoint > reducible;
	reducible.reserve(reducibleSet.size());
	for ( auto m : reducibleSet )
		reducible.push_back( std::move(m.first) );

	//This computes the mappings - we actually don't care about the irreducible data here.
	std::vector< RegularSymmetricGrid::GridPoint > irreducible;
	SymmetryReduction< RegularSymmetricGrid::GridPoint >( symmetry_,
			reducible, irreducible,
			redToIrredMap_ , symRedToIrredMap_,
			irredToRedMap_, symIrredToRedMap_ );

	nIrredPts_ = irreducible.size();

	//The symmetry reduction yields an arbitrary definition of 'irreducible' k vectors.
	//If a constrained is supplied, we need to rotate the mappings to account for that.
	if ( constraintIrreducible.size() != 0 )
	{
		assert( constraintIrreducible.size() % 3 == 0 );

		std::vector<int> RedIndices;
		this->get_list_reducible_lattice_point_indices(constraintIrreducible,RedIndices);

		std::vector<int> constraintIrred = RedIndices;
		for ( auto &ci : constraintIrred)
			ci = redToIrredMap_[ci];

		//previousIrred must be complete and unique or we have a problem
		//a.k.a we cannot fulfill the constraint
		std::set<int> checkIrred(constraintIrred.begin(),constraintIrred.end());
		if ( checkIrred.size() < constraintIrred.size() )
			throw std::runtime_error( "Cannot fulfill irreducible zone constraint:"
					" supplied vectors have overlapping stars" );

		std::set<int> completeIrred(redToIrredMap_.begin(),redToIrredMap_.end());
		if ( checkIrred.size() != completeIrred.size() )
			throw std::runtime_error( "Cannot fulfill irreducible zone constraint:"
					" supplied vectors are not complete in the irreducible zone" );

		//Only if the constraint is not already fulfilled, we need to work ...
		std::set<int> redSet(RedIndices.begin(),RedIndices.end());
		if ( ! (checkIrred == redSet) )
		{
			//walk through the irreducible indices and redirect redToIrredMap_
			//look up which symmetry operation in the star takes one to the constraint grid vector
			for ( int irr = 0 ; irr < nIrredPts_; ++irr)
			{
				int isymNewId = -1;
				for ( int istar = 0 ; istar < symIrredToRedMap_[irr].size(); ++istar)
				{
					auto ret = redSet.find( irredToRedMap_[irr][istar] );
					if ( ret != redSet.end() )
					{
						//irredToRedMap_[irr][istar] is the new irreducible point.
						//Thus apply the inverse to symIrredToRedMap_[irr][istar] to all members of the star.
						isymNewId =  symIrredToRedMap_[irr][istar];
						break;
					}
				}
				assert ( isymNewId >= 0 );
				int invSym = symmetry_.get_index_inverse(isymNewId);
				for ( int istar = 0 ; istar < int(symIrredToRedMap_[irr].size()); ++istar)
				{
					symIrredToRedMap_[irr][istar] = symmetry_.get_group_product(symIrredToRedMap_[irr][istar],invSym);
					int iRed = irredToRedMap_[irr][istar];
					symRedToIrredMap_[iRed] = symmetry_.get_index_inverse(symIrredToRedMap_[irr][istar]);
				}
			}
		}
	}
}

int
RegularSymmetricGrid::get_np_red() const
{
	return this->get_num_points();
}

int
RegularSymmetricGrid::get_np_irred() const
{
	return nIrredPts_;
}

void
RegularSymmetricGrid::get_list_lattice_point_indices(
		std::vector<double> const & points,
		std::vector<int> & irreducibleIndices) const
{
	this->get_list_reducible_lattice_point_indices(points,irreducibleIndices);
	for ( auto &i : irreducibleIndices )
		i = redToIrredMap_[i];
}

void
RegularSymmetricGrid::convert_reducible_irreducible(
		std::vector<int> const& reducibleIndices,
		std::vector<int> & irreducibleIndices) const
{
	if ( irreducibleIndices.size() != reducibleIndices.size() )
		irreducibleIndices.resize(reducibleIndices.size());

	for ( size_t i = 0 ; i < reducibleIndices.size(); ++i )
	{
		assert(  reducibleIndices[i] < this->get_np_red() );
		irreducibleIndices[i] = redToIrredMap_[ reducibleIndices[i] ];
	}
}

std::vector< std::vector<int> > const &
RegularSymmetricGrid::get_maps_irreducible_to_reducible() const
{
	return irredToRedMap_;
}

std::vector< std::vector<int> > const &
RegularSymmetricGrid::get_maps_sym_irred_to_reducible() const
{
	return symIrredToRedMap_;
}

std::vector< int > const &
RegularSymmetricGrid::get_maps_sym_red_to_irreducible() const
{
	return symRedToIrredMap_;
}

std::vector< int > const &
RegularSymmetricGrid::get_maps_red_to_irreducible() const
{
	return redToIrredMap_;
}

bool
RegularSymmetricGrid::is_reci() const
{
	return symmetry_.is_reci();
}

LatticeStructure::Symmetry const &
RegularSymmetricGrid::get_symmetry() const
{
	return symmetry_;
}

} /* namespace LatticeStructure */
} /* namespace elephon */
