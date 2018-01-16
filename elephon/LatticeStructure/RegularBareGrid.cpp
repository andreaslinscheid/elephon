/*	This file RegularBareGrid.cpp is part of elephon.
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
 *  Created on: Jul 7, 2017
 *      Author: A. Linscheid
 */

#include "LatticeStructure/RegularBareGrid.h"
#include "Algorithms/helperfunctions.hpp"
#include <assert.h>
#include <cmath>
#include <algorithm>
#include <set>
#include <stdexcept>
#include <iostream>

namespace elephon
{
namespace LatticeStructure
{

void
RegularBareGrid::initialize(
		std::vector<int> dim,
		bool isReciprocal,
		double gridPrec,
		std::vector<double> shift,
		LatticeModule lattice)
{
	isReciprocalGrid_ = isReciprocal;
	gridPrec_ = gridPrec;
	pointMesh_ = std::move(dim);
	pointShift_ = std::move(shift);
	lattice_ = std::move(lattice);
	assert( pointShift_.size() == 3 );
	assert( pointMesh_.size() == 3 );
	for ( int i = 0 ; i < 3 ; ++i)
	{
		if ( pointShift_[i] >= 1.0 )
			throw std::runtime_error("Not accepting grid shifts by more than one lattice point.");
		pointShift_[i] /= double(pointMesh_[i]);
	}
	numPoints_ = pointMesh_[0]*pointMesh_[1]*pointMesh_[2];
}

int
RegularBareGrid::get_num_points() const
{
	return numPoints_;
}

void
RegularBareGrid::construct_grid_vector_set(
		std::map< RegularBareGrid::GridPoint, int > & reducibleSet) const
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

LatticeStructure::LatticeModule const &
RegularBareGrid::get_lattice() const
{
	return lattice_;
}

int
RegularBareGrid::get_xyz_to_reducible(std::vector<int> const & xyzTouple) const
{
	assert( xyzTouple.size() == 3 );
	assert( (xyzTouple[0] >= 0) && (xyzTouple[0] < pointMesh_[0]) );
	assert( (xyzTouple[1] >= 0) && (xyzTouple[1] < pointMesh_[1]) );
	assert( (xyzTouple[2] >= 0) && (xyzTouple[2] < pointMesh_[2]) );
	return (xyzTouple[2]*pointMesh_[1]+xyzTouple[1])*pointMesh_[0]+xyzTouple[0];
}

int
RegularBareGrid::get_xyz_to_reducible_periodic(std::vector<int> xyzTouple) const
{
	assert( xyzTouple.size() == 3 );
	for ( int i = 0 ; i < 3 ; ++i)
	{
		xyzTouple[i] = xyzTouple[i] < 0 ? xyzTouple[i]+pointMesh_[i] : xyzTouple[i];
		xyzTouple[i] = xyzTouple[i] >= pointMesh_[i] ? xyzTouple[i]-pointMesh_[i] : xyzTouple[i];
	}
	return this->get_xyz_to_reducible(xyzTouple);
}

std::vector<int>
RegularBareGrid::get_reducible_to_xyz(int i) const
{
	std::vector<int> tuple(3);
	this->get_reducible_to_xyz(i, tuple);
	return tuple;
}

void
RegularBareGrid::get_reducible_to_xyz(int i, std::vector<int> & xyz) const
{
	assert(xyz.size() == 3);
	//reverse engineer get_xyz_to_reducible using integer division
	//i = (xyzTouple[2]*pointMesh_[1]+xyzTouple[1])*pointMesh_[0]+xyzTouple[0];
	std::fill(xyz.begin(), xyz.end(), i);
	xyz[2] /= pointMesh_[1]*pointMesh_[0];
	xyz[1] -= xyz[2]*pointMesh_[1]*pointMesh_[0];
	xyz[0] -= xyz[2]*pointMesh_[1]*pointMesh_[0];
	xyz[1] /= pointMesh_[0];
	xyz[0] -= xyz[1]*pointMesh_[0];
	assert( i == xyz[0]+pointMesh_[0]*(xyz[1]+pointMesh_[1]*xyz[2]) );
}

void
RegularBareGrid::get_list_reducible_lattice_point_indices(
		std::vector<double> const & points,
		std::vector<int> & reducibleIndices) const
{
	assert( points.size() % 3 == 0 );

	//Construct the mesh for fast point lookup
	std::map< RegularBareGrid::GridPoint , int > reducibleSet;
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

std::vector<int> const &
RegularBareGrid::get_grid_dim() const
{
	return pointMesh_;
}

std::vector<double>
RegularBareGrid::get_grid_shift() const
{
	std::vector<double> scaled = pointShift_;
	for (int id = 0 ; id < 3;  ++id)
		scaled[id] *= pointMesh_[id];
	return scaled;
}

double
RegularBareGrid::get_grid_prec() const
{
	return gridPrec_;
}

std::vector<double>
RegularBareGrid::get_vector_direct(int i) const
{
	std::vector<int> xyz(3);
	std::vector<double> v(3);
	this->get_vector_direct(i, xyz, v);
	return v;
}

void
RegularBareGrid::get_vector_direct(
		int i,
		std::vector<int> & xyz,
		std::vector<double> & v) const
{
	assert(xyz.size() == 3);
	assert(v.size() == 3);
	this->get_reducible_to_xyz(i, xyz);
	v[0] = double(xyz[0])/double(pointMesh_[0]) + pointShift_[0];
	v[1] = double(xyz[1])/double(pointMesh_[1]) + pointShift_[1];
	v[2] = double(xyz[2])/double(pointMesh_[2]) + pointShift_[2];
	for ( auto &xi : v)
		xi -= std::floor(xi + 0.5);
}

std::vector<double>
RegularBareGrid::get_vectors_direct(
		std::vector<int> const & reducibleIndices) const
{
	const int np = reducibleIndices.size();
	std::vector<double> vectors(np*3);
	std::vector<int> xyz(3);
	std::vector<double> v(3);
	for (int ip = 0 ; ip < np ; ++ip)
	{
		this->get_vector_direct(reducibleIndices[ip], xyz, v);
		for (int i = 0 ; i < 3 ; ++i)
			vectors[ip*3+i] = v[i];
	}
	return vectors;
}

std::vector<double>
RegularBareGrid::get_all_vectors_grid() const
{
	std::vector<int> allGridIndices(numPoints_);
	for (int i = 0 ; i < numPoints_; ++i)
		allGridIndices[i] = i;
	return this->get_vectors_direct(allGridIndices);
}

std::vector<int>
RegularBareGrid::compute_reducible_cube_indices_surrounding_nongrid_point(
		std::vector<double> const& v) const
{
	assert( v.size() == 3 );

	std::vector<int> result(8);
	std::vector<std::pair<int,int> > minMaxEachDim(3);
	for ( int i = 0 ; i < 3; i++ )
	{
		double vfbz = v[i]-pointShift_[i];
		vfbz -= std::floor(vfbz);
		int cellstart = int(std::floor(vfbz*pointMesh_[i]));
		int cellend = cellstart + 1;
		//Apply periodicity: due to finite accuracy of floating point numbers, we have
		// to do this for cellstart, too.
		cellstart = cellstart >= pointMesh_[i] ? cellstart - pointMesh_[i] : cellstart;
		cellend = cellend >= pointMesh_[i] ? cellend - pointMesh_[i] : cellend;
		minMaxEachDim[i] = std::make_pair(cellstart,cellend);
	}

	//3D counter clock wise starting at min min bottom, then clock wise top
	// 1 = min, min, min
	std::vector<int> dima =
		{ minMaxEachDim[0].first, minMaxEachDim[1].first, minMaxEachDim[2].first };
	result[0] = this->get_xyz_to_reducible( dima ) ;

	// 2 = max, min , min
	dima = { minMaxEachDim[0].second, minMaxEachDim[1].first, minMaxEachDim[2].first };
	result[1] = this->get_xyz_to_reducible( dima ) ;

	// 3 = max, max, min
	dima = { minMaxEachDim[0].second, minMaxEachDim[1].second, minMaxEachDim[2].first };
	result[2] = this->get_xyz_to_reducible( dima ) ;

	// 4 = min, max, min
	dima = { minMaxEachDim[0].first, minMaxEachDim[1].second, minMaxEachDim[2].first };
	result[3] = this->get_xyz_to_reducible( dima ) ;

	// 5 = min, min, max
	dima = { minMaxEachDim[0].first, minMaxEachDim[1].first, minMaxEachDim[2].second };
	result[4] = this->get_xyz_to_reducible( dima ) ;

	// 6 = max, min, max
	dima = { minMaxEachDim[0].second, minMaxEachDim[1].first, minMaxEachDim[2].second };
	result[5] = this->get_xyz_to_reducible( dima ) ;

	// 7 = max, max, max
	dima = { minMaxEachDim[0].second, minMaxEachDim[1].second, minMaxEachDim[2].second };
	result[6] = this->get_xyz_to_reducible( dima ) ;

	// 8 = min, max, max
	dima = { minMaxEachDim[0].first, minMaxEachDim[1].second, minMaxEachDim[2].second };
	result[7] = this->get_xyz_to_reducible( dima ) ;

	return result;
}

void
RegularBareGrid::find_closest_reducible_grid_points(
		std::vector<double> const & nonGridPoints,
		std::vector<int> & reducibleGridIndices) const
{
	assert( nonGridPoints.size()%3 == 0 );
	int nptsA = nonGridPoints.size()/3;
	reducibleGridIndices.resize(nptsA);
	for ( int ip = 0 ; ip < nptsA ; ++ip)
	{
		double x = (nonGridPoints[ip*3+0] - std::floor(nonGridPoints[ip*3+0]) - pointShift_[0])*pointMesh_[0];
		double y = (nonGridPoints[ip*3+1] - std::floor(nonGridPoints[ip*3+1]) - pointShift_[1])*pointMesh_[1];
		double z = (nonGridPoints[ip*3+2] - std::floor(nonGridPoints[ip*3+2]) - pointShift_[2])*pointMesh_[2];
		std::vector<int> xyz{ 	static_cast<int>(std::floor(x+0.5)),
								static_cast<int>(std::floor(y+0.5)),
								static_cast<int>(std::floor(z+0.5))};
		reducibleGridIndices[ip] = this->get_xyz_to_reducible_periodic(xyz);
	}
}

void
RegularBareGrid::compute_grid_cubes_surrounding_nongrid_points(
		std::vector<double> const & nonGridPoints,
		std::vector<int> & nonGridPtToCubeMap,
		std::vector<GridCube> & cubes) const
{
	assert( nonGridPoints.size()%3 == 0 );
	int nptsA = nonGridPoints.size()/3;

	//here we compute all the indices of corner points of cubes in the regular grid
	// and accumulate all irregular points in their respective cube
	std::set< GridCube > cubeSet;
	auto lastUsed = cubeSet.end();
	for ( int ip = 0 ; ip < nptsA; ++ip)
	{
		GridCube c( this->compute_reducible_cube_indices_surrounding_nongrid_point(
				std::vector<double>(&nonGridPoints[ip*3],&nonGridPoints[ip*3]+3) ) );
		// we hint the previous location for insert since we expect the irregular grid points
		//	to be close together
		auto ret = cubeSet.insert(lastUsed,c);
		ret->containedIrregularPts_.push_back(ip);
		lastUsed = ret;
	}
	cubes = std::vector< GridCube >( cubeSet.begin(), cubeSet.end() );

	//Now we invert the mapping and fill the vector that maps from nonGridPtToCubeMap
	nonGridPtToCubeMap.resize(nptsA);
	for ( int ic = 0 ; ic < cubes.size(); ++ic)
		for ( auto i : cubes[ic].containedIrregularPts_)
			nonGridPtToCubeMap[i] = ic;
}

std::vector<int>
RegularBareGrid::interpret_fft_dim_input(std::vector<int> fftDim) const
{
	if ( fftDim.size() == 1 )
	{
		int scale = fftDim.at(0);
		fftDim = pointMesh_;
		if ( scale != 0 )
			for ( auto &d : fftDim )
				d *= scale;
	}
	else
	{
		if ( fftDim.size() != 3 )
			throw std::runtime_error("Incorrect grid dimension for fft dimensions.");
		for ( int id = 0 ; id < 3; ++id)
			fftDim[id] = fftDim[id] == 0 ? pointMesh_[id] : fftDim[id];
	}
	return fftDim;
}

void
RegularBareGrid::get_grid_cubes(std::vector<GridCube> & cubes) const
{
	cubes.reserve( pointMesh_[0]*pointMesh_[1]*pointMesh_[2] );
	for ( int iz = 0 ; iz < pointMesh_[2]; ++iz )
		for ( int iy = 0 ; iy < pointMesh_[1]; ++iy )
			for ( int ix = 0 ; ix < pointMesh_[0]; ++ix )
				// consult compute_reducible_cube_indices_surrounding_nongrid_point for the order
				cubes.push_back(
					std::move(GridCube({
				 	 	 	 	 	this->get_xyz_to_reducible_periodic({ix,   iy,   iz}),     // 1 = min, min, min
									this->get_xyz_to_reducible_periodic({ix+1, iy,   iz}),     // 2 = max, min, min
									this->get_xyz_to_reducible_periodic({ix+1, iy+1, iz}),     // 3 = max, max, min
									this->get_xyz_to_reducible_periodic({ix,   iy+1, iz}),     // 4 = min, max, min
									this->get_xyz_to_reducible_periodic({ix,   iy,   iz+1}),   // 5 = min, min, max
									this->get_xyz_to_reducible_periodic({ix+1, iy,   iz+1}),   // 6 = max, min, max
									this->get_xyz_to_reducible_periodic({ix+1, iy+1, iz+1}),   // 7 = max, max, max
									this->get_xyz_to_reducible_periodic({ix,   iy+1, iz+1}),   // 8 = min, max, max
								})));
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
	gridPrec_ = gridPrec;
	this->map_back();
}

void
GridPt::transform( symmetry::SymmetryOperation const & sop )
{
	for (auto &xi : coords_)
		xi -= std::floor(xi+0.5);
	sop.apply(coords_,true);
	this->map_back();
}

std::vector<double>
GridPt::get_coords() const
{
	auto res = coords_;
	for (auto & pi : res)
		pi -= std::floor(pi+0.5);
	return res;
}

bool
operator< (GridPt const& d1, GridPt const& d2)
{
	//Compare position
	assert( d1.gridPrec_ == d2.gridPrec_ );
	for ( int i = 3 ; i--; )
		if ( std::abs( d1.coords_[i] - d2.coords_[i]) > d1.gridPrec_ )
			return d1.coords_[i] < d2.coords_[i];
	return false;
}

void
GridPt::map_back()
{
	for (auto &xi : coords_)
	{
		xi -= std::floor(xi);
		// handle the case where xi is almost equivalent to zero.
		if ( (1.0-xi) < gridPrec_ )
			xi = 0.0;
	}
}

bool
operator< (GridCubeImpl const& d1, GridCubeImpl const& d2)
{
	//Compare position the first corner point. The others must follow the same ordering.
	//This only works for cubes of the same grid size! If cubes d1 and d2 have different size
	//this ordering returns non-sense but there is no way to check this without recording the absolute
	//size of the corner point vectors.
	return d1.cornerIndices_[0] < d2.cornerIndices_[0];
}
} /* namespace detail */

} /* namespace LatticeStructure */
} /* namespace elephon */
