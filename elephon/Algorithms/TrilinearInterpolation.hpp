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
#include <algorithm>

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

	pointsData.resize(gridnum*nB);

	std::vector<double> kcell1(3), kcell2(3);
	int const numTetras = tetraIndexToDataPoints_.size()/4;
	for ( int itetra = 0; itetra < numTetras; ++itetra )
	{
		const int nVThisTetra = tetraContainedIndicesList_[itetra].second.size();

		auto g1it = &( gridDataForRequiredIndices[ nB*tetraIndexToDataPoints_[itetra*4+0]] );
		auto g2it = &( gridDataForRequiredIndices[ nB*tetraIndexToDataPoints_[itetra*4+1]] );
		auto g3it = &( gridDataForRequiredIndices[ nB*tetraIndexToDataPoints_[itetra*4+2]] );
		auto g4it = &( gridDataForRequiredIndices[ nB*tetraIndexToDataPoints_[itetra*4+3]] );

		// construct all the Barycentric coordinates
		std::vector<double> vectorsBarycentric(nVThisTetra*4);
		std::vector<double> vectorsInTetra(nVThisTetra*3);
		std::vector<bool> isInTetra;
		for ( int ik = 0 ; ik < nVThisTetra; ++ik )
		{
			int ikirred = tetraContainedIndicesList_[itetra].second[ik];
			std::copy(&listOfPoints_[ikirred*3], &listOfPoints_[ikirred*3]+3, &vectorsInTetra[ik*3]);
		}
		tetraContainedIndicesList_[itetra].first.check_vectors_inside(
				vectorsInTetra,
				isInTetra,
				vectorsBarycentric);
		assert(std::all_of(isInTetra.begin(), isInTetra.end(), [] (bool a){return a;}));

		for ( int ik = 0 ; ik < nVThisTetra; ++ik )
		{
			int ikirred = tetraContainedIndicesList_[itetra].second[ik];
			for ( int ib= 0 ; ib < nB; ++ib)
			{
				pointsData[ikirred*nB+ib] = vectorsBarycentric[ik*4+0]*g1it[ib] +
											vectorsBarycentric[ik*4+1]*g2it[ib] +
											vectorsBarycentric[ik*4+2]*g3it[ib] +
											vectorsBarycentric[ik*4+3]*g4it[ib] ;
				//check for NaN in debug mode
				assert( pointsData[ikirred*nB+ib] == pointsData[ikirred*nB+ib]);
			}
		}
	}
}

} /* namespace Algorithms */
} /* namespace elephon */


#endif /* ELEPHON_ALGORITHMS_TRILINEARINTERPOLATION_HPP_ */
