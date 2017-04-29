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
#include <complex>
#include <assert.h>

namespace elephon
{
namespace Algorithms
{

namespace detail
{
template<typename T>
struct cmplxTypeTrait
{
	typedef T realT;
};

template<typename T>
struct cmplxTypeTrait<std::complex<T> >
{
	typedef T realT;
};
} /*namespace detail*/

template<typename T>
void TrilinearInterpolation::interpolate(
	size_t nDataPtsPerPoint,
	std::vector<T> const& gridDataForRequiredIndices,
	std::vector<T> & pointsData) const
{
	typedef typename detail::cmplxTypeTrait<T>::realT realT;
	size_t gridnum = listOfPoints_.size()/3;
	assert( gridDataForRequiredIndices.size() == nDataPtsPerPoint*conseqPtsRegularGrid_.size());
	size_t nB = nDataPtsPerPoint;

	typedef std::vector<double> kvect ;
	auto k = kvect(3);
	auto fill_k = [&] (size_t irredKIndex)
	{
		assert( irredKIndex < gridnum );
		auto it = listOfPoints_.begin()+irredKIndex*3;
		std::copy( it, it+3, k.begin() );
	};

	//Grid logic to give the vector made of just the largest and smallest values in every dimension for a cube
	auto get_cell_vectors = [&] ( kvect const& k,kvect & lowerCorner, kvect & upperCorner )
	{
		assert( (upperCorner.size() == lowerCorner.size()) && (lowerCorner.size() == 3) );
		for ( size_t i = 0 ; i < 3; i++ )
		{
			//min, max index in a periodic grid
			double vfbz = k[i]-std::floor(k[i]);
			double vcellmin = std::floor(vfbz*grid_[i])/grid_[i];
			double vcellmax = (std::floor(vfbz*grid_[i])+1.0)/grid_[i];
			lowerCorner[i]=vcellmin;
			upperCorner[i]=vcellmax;
		}
	};

	//Linear interpolation logic
	auto bilinear_interpol = [] (
			double x, double y,
			T f00, T f10, T f11, T f01)
	{
		T a = f00 * realT(1. - x) + f10 * realT(x);
		T b = f01 * realT(1. - x) + f11 * realT(x);
		return a * realT(1. - y) + b * realT(y);
	};

	auto trilinear_interpol = [&] (
			double x, double y, double z,
			T f000, T f100, T f110, T f010,
			T f001, T f101, T f111, T f011)
	{
		auto e = bilinear_interpol(x, y, f000, f100, f110, f010);
		auto f = bilinear_interpol(x, y, f001, f101, f111, f011);
		return e * realT( 1. - z) + f * realT(z);
	};

	if ( pointsData.size() != gridnum*nB )
		pointsData = std::vector<double>(gridnum*nB);

	std::vector<double> kcell1(3), kcell2(3);
	for ( size_t icube = 0; icube < usedGridCubes_.size(); ++icube )
	{
		//use one representative k vector to compute the cell vectors
		fill_k( usedGridCubes_[icube].containedIrregularPts_.front() );
		get_cell_vectors(k,kcell1, kcell2);

		auto g1it = &( gridDataForRequiredIndices[ nB*usedGridCubes_[icube].cornerPoints_[0]] );
		auto g2it = &( gridDataForRequiredIndices[ nB*usedGridCubes_[icube].cornerPoints_[1]] );
		auto g3it = &( gridDataForRequiredIndices[ nB*usedGridCubes_[icube].cornerPoints_[2]] );
		auto g4it = &( gridDataForRequiredIndices[ nB*usedGridCubes_[icube].cornerPoints_[3]] );
		auto g5it = &( gridDataForRequiredIndices[ nB*usedGridCubes_[icube].cornerPoints_[4]] );
		auto g6it = &( gridDataForRequiredIndices[ nB*usedGridCubes_[icube].cornerPoints_[5]] );
		auto g7it = &( gridDataForRequiredIndices[ nB*usedGridCubes_[icube].cornerPoints_[6]] );
		auto g8it = &( gridDataForRequiredIndices[ nB*usedGridCubes_[icube].cornerPoints_[7]] );

		for ( size_t ikirred : usedGridCubes_[icube].containedIrregularPts_ )
		{
			fill_k(ikirred);
			double x = (k[0] - kcell1[0])/(kcell2[0] - kcell1[0]);
			double y = (k[1] - kcell1[1])/(kcell2[1] - kcell1[1]);
			double z = (k[2] - kcell1[2])/(kcell2[2] - kcell1[2]);
			for ( size_t ib= 0 ; ib < nB; ++ib)
			{
				pointsData[ikirred*nB+ib] = trilinear_interpol(
						x , y, z,
						g1it[ib], g2it[ib],
						g3it[ib], g4it[ib],
						g5it[ib], g6it[ib],
						g7it[ib], g8it[ib] );
				//check for NaN in debug mode
				assert( pointsData[ikirred*nB+ib] == pointsData[ikirred*nB+ib]);
			}
		}
	}
}

} /* namespace Algorithms */
} /* namespace elephon */


#endif /* ELEPHON_ALGORITHMS_TRILINEARINTERPOLATION_HPP_ */
