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
		std::shared_ptr<const LatticeStructure::TetrahedraGrid> grid)
				: tetraGrid_(grid)
{

}

void
TrilinearInterpolation::data_query(
		std::vector<double> listOfPoints,
		std::vector<int> & requiredGridIndices)
{
	listOfPoints_ = std::move(listOfPoints);
	//Map back to the periodic cell [0,1[ this is done, because the tetrahedra are defined in this cell.
	for ( auto & xi : listOfPoints_ )
		xi -= std::floor(xi);
	int nPts = listOfPoints_.size()/3;

	std::map<LatticeStructure::Tetrahedron, std::vector<int>> tetraContainedIndicesListMap;
	tetraGrid_->compute_grid_tetrahedra_surrounding_nongrid_points(
			listOfPoints_,
			tetraContainedIndicesListMap);

	tetraContainedIndicesList_.clear();
	tetraContainedIndicesList_.reserve(tetraContainedIndicesListMap.size());
	for ( auto const &tm : tetraContainedIndicesListMap)
		tetraContainedIndicesList_.push_back(tm);

	//Now we determine all appearing regular grid indices
	std::set<int> conseqRegularGridIndices;
	for ( auto tetra : tetraContainedIndicesList_)
		conseqRegularGridIndices.insert(tetra.first.get_corner_indices().begin(),
										tetra.first.get_corner_indices().end() );

	requiredGridIndices =  std::vector<int>(conseqRegularGridIndices.begin(),
											conseqRegularGridIndices.end());
	conseqRegularGridIndices.clear();

	//If data interpolation is called, we need
	//	1) to associate with each data index a given regular grid index
	//	2) store for each cell index, the data index of each corner point

	//	1) Enumerate points (data layout) and map each data index 'data_i' to a regular grid index
	conseqPtsRegularGrid_.resize( requiredGridIndices.size() );
	//The following variable is used to reset the regular grid index in the cubes to the data layout
	std::map<int,int> invConseqPtsRegularGrid;
	int i = 0;
	for ( auto ig : requiredGridIndices )
	{
		invConseqPtsRegularGrid[ig] = i; // guaranteed to be unique by map
		conseqPtsRegularGrid_[i++]= ig;
	}

	//	2) Set for each tetrahedron in tetraContainedIndicesList_ the corner indices to the ones
	//		in the data array. The result is stored in  tetraIndexToDataPoints_.
	// 		Also cross check that we have mapped all points.
	tetraIndexToDataPoints_.resize(4*tetraContainedIndicesList_.size());
	int nPtsCheck = 0;
	int itetra = 0;
	for ( auto & tetra : tetraContainedIndicesList_)
	{
		std::copy(tetra.first.get_corner_indices().begin(),
				  tetra.first.get_corner_indices().end(),
				  &tetraIndexToDataPoints_[itetra*4]);
		for ( auto i : tetra.second)
		{
			auto ret = conseqRegularGridIndices.insert(i);
			if( not ret.second )
				throw std::logic_error(std::string("Index of the non-grid points #")
					+std::to_string(i)+" appears multiple times. Internal logic error!");
		}
		nPtsCheck += tetra.second.size();
		for (int i = 0 ; i < 4 ; ++i)
		{
			auto it = invConseqPtsRegularGrid.find( tetraIndexToDataPoints_[itetra*4+i] );
			if ( it == invConseqPtsRegularGrid.end() )
				throw std::logic_error ("Internal programming error, !");
			tetraIndexToDataPoints_[itetra*4+i] = it->second;
		}
		itetra++;
	}

	if ( nPts != nPtsCheck)
		throw std::logic_error ("We have lost or gained some grid points. Internal programming error!");
}

void TrilinearInterpolation::interpolate(
		int nDataPtsPerPoint,
		std::vector<double> const& gridDataForRequiredIndices,
		std::vector<double> & pointsData) const
{
	this->interpolate<double>(nDataPtsPerPoint,gridDataForRequiredIndices,pointsData);
}

} /* namespace Algorithms */
} /* namespace elephon */
