/*	This file TrilinearInterpolation.cpp is part of elephon.
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
 *  Created on: Apr 28, 2017
 *      Author: A. Linscheid
 */

#include "TrilinearInterpolation.h"
#include <assert.h>
#include <set>
#include <exception>
#include <cmath>
#include <map>

namespace elephon
{
namespace Algorithms
{

TrilinearInterpolation::TrilinearInterpolation(
		std::vector<size_t>  grid) : grid_(std::move(grid))
{
	assert( grid_.size() == 3 );
}

void TrilinearInterpolation::data_query(
		std::vector<double> listOfPoints,
		std::vector<size_t> & requiredGridIndices)
{
	listOfPoints_ = std::move(listOfPoints);
	//Map back to the periodic cell!
	for ( auto & xi : listOfPoints_ )
		xi -= std::floor(xi);
	size_t nPts = listOfPoints_.size()/3;

	//here we compute all the indices of corner points of cubes in the regular grid
	// and accumulate all irregular grid points in their respective cube
	// and save the cubes for later use
	std::set< GridCube > cubeSet;
	std::vector<double> v(3);
	auto lastUsed = cubeSet.end();
	for ( size_t ip = 0 ; ip < nPts ; ++ip)
	{
		auto it = listOfPoints_.begin()+ip*3;
		std::copy(it,it+3,v.begin() );
		GridCube c( this->get_cube_indices_surrounding(v) );
		// we hint the previous location for insert since we expect the irregular grid points
		//	to be close together
		auto ret = cubeSet.insert(lastUsed,c);
		ret->containedIrregularPts_.push_back(ip);
		lastUsed = ret;
	}
	usedGridCubes_ = std::vector< GridCube >( cubeSet.begin(), cubeSet.end() );

	//Now we determine all appearing regular grid indices
	std::set<size_t> conseqRegularGridIndices;
	for ( auto cube : usedGridCubes_)
		conseqRegularGridIndices.insert( cube.cornerPoints_.begin(), cube.cornerPoints_.end() );

	requiredGridIndices =  std::vector<size_t>(
			conseqRegularGridIndices.begin(),
			conseqRegularGridIndices.end());

	//If data interpolation is called, we need
	//	1) to associate with each data index a given regular grid index
	//	2) store for each cell index, the data index of each corner point

	//	1) Enumerate points (data layout) and map each data index 'data_i' to a regular grid index
	conseqPtsRegularGrid_ = std::vector<size_t>( conseqRegularGridIndices.size() );
	//The following variable is used to reset the regular grid index in the cubes to the data layout
	std::map<size_t,size_t> invConseqPtsRegularGrid;
	size_t i = 0;
	for ( auto ig : conseqRegularGridIndices )
	{
		invConseqPtsRegularGrid[ig] = i; // guaranteed to be unique by set
		conseqPtsRegularGrid_[i++]= ig;
	}

	//	2) Reset for each cube in usedGridCubes_ the corner index to the index
	//		in the data array. Also cross check that we have mapped all points.
	size_t nPtsCheck = 0;
	for ( auto & cube : usedGridCubes_)
	{
		for ( auto & coi : cube.cornerPoints_ )
		{
			auto it = invConseqPtsRegularGrid.find( coi );
			if ( it == invConseqPtsRegularGrid.end() )
				throw std::logic_error ("Internal programming error, !");
			coi = it->second;
		}
		nPtsCheck += cube.containedIrregularPts_.size();
	}

	if ( nPts != nPtsCheck)
		throw std::logic_error ("We have lost some grid points. Internal programming error!");
}

void TrilinearInterpolation::interpolate(
		size_t nDataPtsPerPoint,
		std::vector<double> const& gridDataForRequiredIndices,
		std::vector<double> & pointsData) const
{
	this->interpolate<double>(nDataPtsPerPoint,gridDataForRequiredIndices,pointsData);
}

bool TrilinearInterpolation::GridCube::operator< (GridCube const& other) const
{
	for ( size_t i = 0 ; i < cornerPoints_.size(); ++i)
	{
		if ( cornerPoints_[i] < other.cornerPoints_[i] )
			return true;
		if ( other.cornerPoints_[i] < cornerPoints_[i] )
			return false;
	}
	return false;
};

std::vector<size_t>
TrilinearInterpolation::get_cube_indices_surrounding(
		std::vector<double> const& v) const
{
	assert( v.size()== 3 );

	std::vector<size_t> result(8,0);
	std::vector<std::pair<size_t,size_t> > minMaxEachDim;
	for ( size_t i = 0 ; i < 3; i++ )
	{
		//min, max index in a grid_ periodic grid
		double vfbz = v[i]-std::floor(v[i]);
		int cellstart = static_cast<int>(std::floor(vfbz*grid_[i]));
		int cellend = cellstart + 1;
		cellstart -= (cellstart/grid_[i])*grid_[i];
		cellend -= (cellend/grid_[i])*grid_[i];
		minMaxEachDim.push_back(
				std::make_pair(static_cast<size_t>(cellstart),static_cast<size_t>(cellend)));
	}

	auto xyz_to_conseq = [&] (std::vector<size_t> const& tuple )
		{
			auto conseq_index = tuple.front();
			for ( int i = 1; i < 3 ; ++i)
			{
				conseq_index *= grid_[i];
				conseq_index += tuple[i];
			}
			return conseq_index;
		};

	//3D counter clock wise starting at min min bottom, then clock wise top
	// 1 = min, min, min
	std::vector<size_t> dima =
		{ minMaxEachDim[0].first, minMaxEachDim[1].first, minMaxEachDim[2].first };
	result[0] =  xyz_to_conseq( dima ) ;

	// 2 = max, min , min
	dima = { minMaxEachDim[0].second, minMaxEachDim[1].first, minMaxEachDim[2].first };
	result[1] = xyz_to_conseq( dima ) ;

	// 3 = max, max, min
	dima = { minMaxEachDim[0].second, minMaxEachDim[1].second, minMaxEachDim[2].first };
	result[2] = xyz_to_conseq( dima );

	// 4 = min, max, min
	dima = { minMaxEachDim[0].first, minMaxEachDim[1].second, minMaxEachDim[2].first };
	result[3] = xyz_to_conseq( dima );

	// 5 = min, min, max
	dima = { minMaxEachDim[0].first, minMaxEachDim[1].first, minMaxEachDim[2].second };
	result[4] = xyz_to_conseq( dima );

	// 6 = max, min, max
	dima = { minMaxEachDim[0].second, minMaxEachDim[1].first, minMaxEachDim[2].second };
	result[5] = xyz_to_conseq( dima );

	// 7 = max, max, max
	dima = { minMaxEachDim[0].second, minMaxEachDim[1].second, minMaxEachDim[2].second };
	result[6] = xyz_to_conseq( dima );

	// 8 = min, max, max
	dima = { minMaxEachDim[0].first, minMaxEachDim[1].second, minMaxEachDim[2].second };
	result[7] = xyz_to_conseq( dima );

	return result;
}

} /* namespace Algorithms */
} /* namespace elephon */
