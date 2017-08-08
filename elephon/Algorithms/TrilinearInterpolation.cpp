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
		LatticeStructure::RegularBareGrid grid)
				: grid_(std::move(grid))
{

}

void TrilinearInterpolation::data_query(
		std::vector<double> listOfPoints,
		std::vector<int> & requiredGridIndices)
{
	listOfPoints_ = std::move(listOfPoints);
	//Map back to the periodic cell!
	for ( auto & xi : listOfPoints_ )
		xi -= std::floor(xi);
	int nPts = listOfPoints_.size()/3;

	std::vector<int> dummy;
	grid_.compute_grid_cubes_surrounding_nongrid_points(
			listOfPoints_,dummy,usedGridCubes_);

	//Now we determine all appearing regular grid indices
	std::set<int> conseqRegularGridIndices;
	for ( auto cube : usedGridCubes_)
		conseqRegularGridIndices.insert( cube.cornerIndices_.begin(), cube.cornerIndices_.end() );

	requiredGridIndices =  std::vector<int>(
			conseqRegularGridIndices.begin(),
			conseqRegularGridIndices.end());

	//If data interpolation is called, we need
	//	1) to associate with each data index a given regular grid index
	//	2) store for each cell index, the data index of each corner point

	//	1) Enumerate points (data layout) and map each data index 'data_i' to a regular grid index
	conseqPtsRegularGrid_.resize( conseqRegularGridIndices.size() );
	//The following variable is used to reset the regular grid index in the cubes to the data layout
	std::map<int,int> invConseqPtsRegularGrid;
	int i = 0;
	for ( auto ig : requiredGridIndices )
	{
		invConseqPtsRegularGrid[ig] = i; // guaranteed to be unique by set
		conseqPtsRegularGrid_[i++]= ig;
	}

	//	2) Reset for each cube in usedGridCubes_ the corner index to the index
	//		in the data array. Also cross check that we have mapped all points.
	int nPtsCheck = 0;
	for ( auto & cube : usedGridCubes_)
	{
		for ( auto & coi : cube.cornerIndices_ )
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

void
TrilinearInterpolation::get_cell_vectors( std::vector<double> const& k,
		std::vector<double> & lowerCorner,
		std::vector<double> & upperCorner ) const
{
	//Grid logic to give the vector made of just the largest and smallest values in every dimension for a cube
	assert( (upperCorner.size() == lowerCorner.size()) && (lowerCorner.size() == 3) );
	for ( size_t i = 0 ; i < 3; i++ )
	{
		double vfbz = k[i]-std::floor(k[i]);
		double vcellmin = std::floor(vfbz*grid_.get_grid_dim()[i])/grid_.get_grid_dim()[i];
		double vcellmax = (std::floor(vfbz*grid_.get_grid_dim()[i])+1.0)/grid_.get_grid_dim()[i];
		lowerCorner[i]=vcellmin;
		upperCorner[i]=vcellmax;
	}
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
