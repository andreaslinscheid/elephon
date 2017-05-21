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
		LatticeModule lattice)
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

	std::set< RegularGrid::GridPoint > reducibleSet;
	this->construct_grid_vector_set( reducibleSet );
	std::vector< RegularGrid::GridPoint > reducible( reducibleSet.begin()  , reducibleSet.end() );

	//This computes the mappings - we actually don't care about the irreducible data here.
	std::vector< RegularGrid::GridPoint > irreducible;
	SymmetryReduction< RegularGrid::GridPoint >( symmetry_,
			reducible, irreducible,
			redToIrredMap_ , symRedToIrredMap_,
			irredToRredMap_, symIrredToRedMap_ );

	nRedPts_ = reducible.size();
	nIrredPts_ = irreducible.size();
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
RegularGrid::get_list_lattice_point_indices(
		std::vector<double> const & points,
		std::vector<int> & irreducibleIndices) const
{
	assert( points.size() % 3 == 0 );

	auto vect_to_touple = [&] (double x, double y, double z) {
		//map to the range [0,1[
		x -= std::floor(x);
		y -= std::floor(y);
		z -= std::floor(z);
		return std::vector<int>( { 	static_cast<int>(std::floor((x-pointShift_[0])*pointMesh_[0]+0.5)),
									static_cast<int>(std::floor((y-pointShift_[1])*pointMesh_[1]+0.5)),
									static_cast<int>(std::floor((z-pointShift_[2])*pointMesh_[2]+0.5))	} );
	};

	//Construct the mesh for fast point lookup
	std::set< RegularGrid::GridPoint > reducibleSet;
	this->construct_grid_vector_set( reducibleSet );

	int nK = static_cast<int>(points.size())/3;
	if ( static_cast<int>(irreducibleIndices.size()) != nK)
		irreducibleIndices.resize( nK );
	for (int ik = 0 ; ik < nK ; ++ik)
	{
		GridPoint gp( std::vector<double>({points[ik*3],points[ik*3+1],points[ik*3+2]}) , gridPrec_ );
		auto it = reducibleSet.find( gp );
		if ( it == reducibleSet.end() )
			throw std::runtime_error( "Grid vectors passed are not matching any grid vector " );

		irreducibleIndices[ik] = redToIrredMap_[ this->get_xyz_to_reducible(
						vect_to_touple(points[ik*3],points[ik*3+1],points[ik*3+2]) ) ];
	}
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

void
RegularGrid::convert_reducible_irreducible(
		std::vector<int> const& reducibleIndices,
		std::vector<int> & irreducibleIndices,
		std::vector<int> & symmetryIndexRedToIrred) const
{
	if ( irreducibleIndices.size() != reducibleIndices.size() )
		irreducibleIndices.resize(reducibleIndices.size());

	if ( symmetryIndexRedToIrred.size() != reducibleIndices.size() )
		symmetryIndexRedToIrred.resize(reducibleIndices.size());

	for ( size_t i = 0 ; i < reducibleIndices.size(); ++i )
	{
		assert(  reducibleIndices[i] < nRedPts_ );
		irreducibleIndices[i] = redToIrredMap_[ reducibleIndices[i] ];
		symmetryIndexRedToIrred[i] = symRedToIrredMap_[ reducibleIndices[i] ];
	}
}

void
RegularGrid::construct_grid_vector_set(
		std::set< RegularGrid::GridPoint > & reducibleSet) const
{
	reducibleSet.clear();
	double d[3] = {1/double(pointMesh_[0]),1/double(pointMesh_[1]),1/double(pointMesh_[2])};
	for (int k = 0 ; k < pointMesh_[2]; ++k)
		for (int j = 0 ; j < pointMesh_[1]; ++j)
			for (int i = 0 ; i < pointMesh_[0]; ++i)
			{
				GridPoint gp( std::vector<double>(
						{ i*d[0]+pointShift_[0],
							j*d[1]+pointShift_[1],
							k*d[2]+pointShift_[2] } ), gridPrec_ );
				auto ret = reducibleSet.insert( std::move(gp) );
				if ( ! ret.second )
					throw std::runtime_error( "Error generating regular grid mesh." );
			}
}

int
RegularGrid::get_xyz_to_reducible(std::vector<int> const & xyzTouple) const
{
	assert( xyzTouple.size() == 3 );
	return (xyzTouple[2]*pointMesh_[1]+xyzTouple[1])*pointMesh_[0]+xyzTouple[0];
}

bool
RegularGrid::is_reci() const
{
	return symmetry_.is_reci();
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
	assert( d1.gridPrec_ == d2.gridPrec_ );

	auto cmp = [&] (double a, double b ) {
		if ( std::fabs(a-b) < d1.gridPrec_ )
			return 0;
		return (a < b ? 1 : -1);
	};

	//Compare position
	for ( int i = 0 ; i < 3; ++i)
		if ( cmp( d1.coords_[i],d2.coords_[i]) )
			return d1.coords_[i] < d2.coords_[i];

	//Position is (almost) equal thus d1<d2 is false
	return false;
}
} /* namespace detail */

} /* namespace LatticeStructure */
} /* namespace elephon */
