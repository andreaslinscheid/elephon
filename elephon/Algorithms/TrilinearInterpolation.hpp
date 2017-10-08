/*	This file TrilinearInterpolation.hpp is part of elephon.
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
 *  Created on: Apr 29, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_ALGORITHMS_TRILINEARINTERPOLATION_HPP_
#define ELEPHON_ALGORITHMS_TRILINEARINTERPOLATION_HPP_

#include "Algorithms/TrilinearInterpolation.h"
#include "Algorithms/helperfunctions.hpp"
#include <complex>
#include <assert.h>
#include <iostream>

namespace elephon
{
namespace Algorithms
{

template<typename T>
void TrilinearInterpolation::interpolate(
	int nDataPtsPerPoint,
	std::vector<T> const& gridDataForRequiredIndices,
	std::vector<T> & pointsData) const
{
	size_t gridnum = listOfPoints_.size()/3;
	assert( gridDataForRequiredIndices.size() == nDataPtsPerPoint*conseqPtsRegularGrid_.size());
	int nB = nDataPtsPerPoint;

	typedef std::vector<double> kvect ;
	auto k = kvect(3);
	auto fill_k = [&] (size_t irredKIndex)
	{
		assert( irredKIndex < gridnum );
		auto it = listOfPoints_.begin()+irredKIndex*3;
		std::copy( it, it+3, k.begin() );
	};

	if ( pointsData.size() != gridnum*nB )
		pointsData = std::vector<double>(gridnum*nB);

	std::vector<double> kcell1(3), kcell2(3);
	for ( int icube = 0; icube < usedGridCubes_.size(); ++icube )
	{
		//use one representative k vector to compute the cell vectors
		fill_k( usedGridCubes_[icube].containedIrregularPts_.front() );
		this->get_cell_vectors(k, kcell1, kcell2);

		auto g1it = &( gridDataForRequiredIndices[ nB*usedGridCubes_[icube].cornerIndices_[0]] );
		auto g2it = &( gridDataForRequiredIndices[ nB*usedGridCubes_[icube].cornerIndices_[1]] );
		auto g3it = &( gridDataForRequiredIndices[ nB*usedGridCubes_[icube].cornerIndices_[2]] );
		auto g4it = &( gridDataForRequiredIndices[ nB*usedGridCubes_[icube].cornerIndices_[3]] );
		auto g5it = &( gridDataForRequiredIndices[ nB*usedGridCubes_[icube].cornerIndices_[4]] );
		auto g6it = &( gridDataForRequiredIndices[ nB*usedGridCubes_[icube].cornerIndices_[5]] );
		auto g7it = &( gridDataForRequiredIndices[ nB*usedGridCubes_[icube].cornerIndices_[6]] );
		auto g8it = &( gridDataForRequiredIndices[ nB*usedGridCubes_[icube].cornerIndices_[7]] );

		for ( int ikirred : usedGridCubes_[icube].containedIrregularPts_ )
		{
			fill_k(ikirred);
			double x = (k[0] - kcell1[0])/(kcell2[0] - kcell1[0]);
			double y = (k[1] - kcell1[1])/(kcell2[1] - kcell1[1]);
			double z = (k[2] - kcell1[2])/(kcell2[2] - kcell1[2]);
			for ( int ib= 0 ; ib < nB; ++ib)
			{
				pointsData[ikirred*nB+ib] =
					helperfunctions::interpolate_single_cube_realtive(
						x , y, z,
						g1it[ib], g2it[ib], g3it[ib], g4it[ib],
						g5it[ib], g6it[ib], g7it[ib], g8it[ib] );
				//check for NaN in debug mode
				assert( pointsData[ikirred*nB+ib] == pointsData[ikirred*nB+ib]);
			}
		}
	}
}

template<typename T>
void
TrilinearInterpolation::interpolate_within_single_cube(
		std::vector<double> const & ptsInCube,
		std::vector<std::vector<T>> const & cornerData,
		std::vector<T> & interpolData) const
{
	assert(ptsInCube.size()%3 == 0);
	assert( cornerData.size() == 8 );
	int nD = cornerData[0].size();
	assert( (nD == cornerData[1].size()) && (nD == cornerData[2].size()) &&
			(nD == cornerData[3].size()) && (nD == cornerData[4].size()) &&
			(nD == cornerData[5].size()) && (nD == cornerData[6].size()) && (nD == cornerData[7].size()) );

	interpolData.resize(nD);
	for (int ip = 0 ; ip < ptsInCube.size()/3 ; ++ip)
	{
		double x = ptsInCube[ip*3+0];
		double y = ptsInCube[ip*3+1];
		double z = ptsInCube[ip*3+2];
		assert( (x >= 0.0) && (x < 1.0) );
		assert( (y >= 0.0) && (y < 1.0) );
		assert( (z >= 0.0) && (z < 1.0) );
		for ( int id = 0 ; id < nD ; ++id )
			interpolData[id] = helperfunctions::interpolate_single_cube_realtive(
					x,y,z,
					cornerData[0][id], cornerData[1][id], cornerData[2][id], cornerData[3][id],
					cornerData[4][id], cornerData[5][id], cornerData[6][id], cornerData[7][id] );
	}
}

} /* namespace Algorithms */
} /* namespace elephon */


#endif /* ELEPHON_ALGORITHMS_TRILINEARINTERPOLATION_HPP_ */
