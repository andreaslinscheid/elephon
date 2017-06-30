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

#include "RegularGrid.h"
#include "SymmetryReduction.h"
#include <assert.h>
#include <cmath>
#include <stdexcept>
#include <iostream>

namespace elephon
{
namespace LatticeStructure
{

void
RegularGrid::initialize(
		double gridPrec,
		std::vector<int> dim,
		std::vector<double> shift,
		Symmetry symmetry,
		LatticeModule lattice,
		std::vector<double> const & constraintIrreducible)
{
	gridPrec_ = gridPrec;
	pointMesh_ = std::move( dim );
	assert( pointMesh_.size() == 3 );
	pointShift_ = std::move( shift );
	assert( pointShift_.size() == 3 );
	for ( int i = 0 ; i < 3 ; ++i)
	{
		if ( pointShift_[i] >= 1.0 )
			throw std::runtime_error("Not accepting grid shifts by more than one lattice point.");
		pointShift_[i] /= double(pointMesh_[i]);
	}
	symmetry_ = std::move( symmetry );
	lattice_ = std::move( lattice );

	std::map< RegularGrid::GridPoint, int > reducibleSet;
	this->construct_grid_vector_set( reducibleSet );
	std::vector< RegularGrid::GridPoint > reducible;
	reducible.reserve(reducibleSet.size());
	for ( auto m : reducibleSet )
		reducible.push_back( std::move(m.first) );

	//This computes the mappings - we actually don't care about the irreducible data here.
	std::vector< RegularGrid::GridPoint > irreducible;
	SymmetryReduction< RegularGrid::GridPoint >( symmetry_,
			reducible, irreducible,
			redToIrredMap_ , symRedToIrredMap_,
			irredToRedMap_, symIrredToRedMap_ );

	nRedPts_ = reducible.size();
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
			throw std::runtime_error( "Cannot fulfill irreducible zone constraint: supplied vectors have overlapping stars" );

		std::set<int> completeIrred(redToIrredMap_.begin(),redToIrredMap_.end());
		if ( checkIrred.size() != completeIrred.size() )
			throw std::runtime_error( "Cannot fulfill irreducible zone constraint:"
					" supplied vectors are not complete in the irreducible zone" );

		//Only if the constraint is not already fulfilled, we need to work ...
		std::set<int> redSet(RedIndices.begin(),RedIndices.end());
		if ( ! (checkIrred == redSet) )
		{
			//walk through the irreducible indices and redicret redToIrredMap_
			//look up which symmetry operation in the star takes one to the constraint grid vector
			for ( int irr = 0 ; irr < nIrredPts_; ++irr)
			{
				int isymNewId = -1;
				for ( int istar = 0 ; istar < int(symIrredToRedMap_[irr].size()); ++istar)
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
RegularGrid::get_np_red() const
{
	return nRedPts_;
}

int
RegularGrid::get_np_irred() const
{
	return nIrredPts_;
}

void
RegularGrid::get_list_reducible_lattice_point_indices(
		std::vector<double> const & points,
		std::vector<int> & reducibleIndices) const
{
	assert( points.size() % 3 == 0 );

	//Construct the mesh for fast point lookup
	std::map< RegularGrid::GridPoint , int > reducibleSet;
	this->construct_grid_vector_set( reducibleSet );

	int nK = points.size()/3;
	reducibleIndices.resize( nK );
	for (int ik = 0 ; ik < nK ; ++ik)
	{
		GridPoint gp( std::vector<double>({points[ik*3],points[ik*3+1],points[ik*3+2]}) , gridPrec_ );
		auto it = reducibleSet.find( gp );
		if ( it == reducibleSet.end() )
			throw std::runtime_error( "Grid vectors passed are not matching any grid vector " );

		reducibleIndices[ik] =  it->second ;
	}
}

void
RegularGrid::get_list_lattice_point_indices(
		std::vector<double> const & points,
		std::vector<int> & irreducibleIndices) const
{
	this->get_list_reducible_lattice_point_indices(points,irreducibleIndices);
	for ( auto &i : irreducibleIndices )
		i = redToIrredMap_[i];
}

void
RegularGrid::convert_reducible_irreducible(
		std::vector<int> const& reducibleIndices,
		std::vector<int> & irreducibleIndices) const
{
	if ( irreducibleIndices.size() != reducibleIndices.size() )
		irreducibleIndices.resize(reducibleIndices.size());

	for ( size_t i = 0 ; i < reducibleIndices.size(); ++i )
	{
		assert(  reducibleIndices[i] < nRedPts_ );
		irreducibleIndices[i] = redToIrredMap_[ reducibleIndices[i] ];
	}
}

std::vector< std::vector<int> > const &
RegularGrid::get_maps_irreducible_to_reducible() const
{
	return irredToRedMap_;
}

std::vector< std::vector<int> > const &
RegularGrid::get_maps_sym_irred_to_reducible() const
{
	return symIrredToRedMap_;
}

std::vector< int > const &
RegularGrid::get_maps_sym_red_to_irreducible() const
{
	return symRedToIrredMap_;
}

std::vector< int > const &
RegularGrid::get_maps_red_to_irreducible() const
{
	return redToIrredMap_;
}

void
RegularGrid::construct_grid_vector_set(
		std::map< RegularGrid::GridPoint, int > & reducibleSet) const
{
	reducibleSet.clear();
	auto hint = reducibleSet.end();
	for ( int i = 0 ; i < pointMesh_[0]*pointMesh_[1]*pointMesh_[2]; ++i)
	{
		GridPoint gp( this->get_vector_direct(i), gridPrec_ );
		hint = reducibleSet.insert(hint, std::move(std::make_pair( std::move(gp), i ) ) );
		if ( (++hint) != reducibleSet.end()  )
			throw std::runtime_error( "Error generating regular grid mesh.\n"
					"New grid point not placed at the end, indicating that ordering "
					"conventions differ which they must not" );
	}
}

int
RegularGrid::get_xyz_to_reducible(std::vector<int> const & xyzTouple) const
{
	assert( xyzTouple.size() == 3 );
	return (xyzTouple[2]*pointMesh_[1]+xyzTouple[1])*pointMesh_[0]+xyzTouple[0];
}

std::vector<int>
RegularGrid::get_reducible_to_xyz(int i) const
{
	//reverse engineer get_xyz_to_reducible using integer division
	//i = (xyzTouple[2]*pointMesh_[1]+xyzTouple[1])*pointMesh_[0]+xyzTouple[0];
	std::vector<int> tuple(3,i);
	tuple[2] /= pointMesh_[1]*pointMesh_[0];
	tuple[1] -= tuple[2]*pointMesh_[1]*pointMesh_[0];
	tuple[0] -= tuple[2]*pointMesh_[1]*pointMesh_[0];
	tuple[1] /= pointMesh_[0];
	tuple[0] -= tuple[1]*pointMesh_[0];
	return tuple;
}

bool
RegularGrid::is_reci() const
{
	return symmetry_.is_reci();
}

LatticeStructure::Symmetry const &
RegularGrid::get_symmetry() const
{
	return symmetry_;
}

std::vector<int> const &
RegularGrid::get_grid_dim() const
{
	return pointMesh_;
}

std::vector<double>
RegularGrid::get_grid_shift() const
{
	std::vector<double> scaled = pointShift_;
	for (int id = 0 ; id < 3;  ++id)
		scaled[id] *= pointMesh_[id];
	return scaled;
}

double
RegularGrid::get_grid_prec() const
{
	return gridPrec_;
}

std::vector<double>
RegularGrid::get_vector_direct(int i) const
{
	auto xyz = this->get_reducible_to_xyz(i);
	std::vector<double> v({ double(xyz[0])/double(pointMesh_[0]) + pointShift_[0],
							double(xyz[1])/double(pointMesh_[1]) + pointShift_[1],
							double(xyz[2])/double(pointMesh_[2]) + pointShift_[2] } );
	for ( auto &xi : v)
		xi -= std::floor(xi + 0.5);
	return v;
}

namespace detail
{

GridPt::GridPt()
{

}

GridPt::GridPt(std::vector<double> coords, double gridPrec)
{
	this->initialize( std::move(coords),gridPrec );
}

void
GridPt::initialize(std::vector<double> coords, double gridPrec)
{
	assert( coords.size() == 3 );
	coords_ = std::move(coords);
	for ( auto &xi : coords_)
		xi -= std::floor( xi + 0.5);
	gridPrec_ = gridPrec;
}

void
GridPt::transform( LatticeStructure::Symmetry::SymmetryOperation const & op)
{
	op.apply( coords_ , /* lattice periodic = */ true );
}

bool
operator< (GridPt const& d1, GridPt const& d2)
{
	//Compare position
	assert( d1.gridPrec_ == d2.gridPrec_ );
	for ( int i = 3 ; i--; )
		if ( std::abs( d1.coords_[i] - d2.coords_[i]) > d1.gridPrec_ )
		{
			double range1 = (d1.coords_[i] < 0 ? d1.coords_[i] + 1.0 : d1.coords_[i]);
			double range2 = (d2.coords_[i] < 0 ? d2.coords_[i] + 1.0 : d2.coords_[i]);
			return range1 < range2;
		}
	return false;
}
} /* namespace detail */

} /* namespace LatticeStructure */
} /* namespace elephon */
