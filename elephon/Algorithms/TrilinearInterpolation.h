/*	This file TrilinearInterpolation.h is part of elephon.
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

#ifndef ELEPHON_ALGORITHMS_TRILINEARINTERPOLATION_H_
#define ELEPHON_ALGORITHMS_TRILINEARINTERPOLATION_H_

#include <vector>
#include <cstdlib>

namespace elephon
{
namespace Algorithms
{

class TrilinearInterpolation
{
public:
	TrilinearInterpolation(
			std::vector<size_t> grid);

	/**
	 * requiredGridIndices will be given in a consecutively ordered grid index
	 * where where x is the slowest running: (ix,iy,iz) <==> (ix*grid[y]+iy)*grid[z]+iz
	 */
	void data_query(
			std::vector<double> listOfPoints,
			std::vector<size_t> & requiredGridIndices);
	/**
	 *
	 */
	void interpolate(
			size_t nDataPtsPerPoint,
			std::vector<double> const& gridDataForRequiredIndices,
			std::vector<double> & pointsData) const;

private:

	struct GridCube
	{
		GridCube( std::vector<size_t> cornerPoints )
				: cornerPoints_( std::move(cornerPoints) ){};

		std::vector<size_t> cornerPoints_;
		mutable std::vector<size_t> containedIrregularPts_;

		bool operator< (GridCube const& other) const;
	};

	std::vector<size_t> grid_;

	std::vector<double> listOfPoints_;

	std::vector<GridCube> usedGridCubes_;

	std::vector<size_t> conseqPtsRegularGrid_;

	std::vector<size_t> get_cube_indices_surrounding(std::vector<double> const& v) const;

	template<typename T>
	void interpolate(
			size_t nDataPtsPerPoint,
			std::vector<T> const& gridDataForRequiredIndices,
			std::vector<T> & pointsData) const;

};

} /* namespace Algorithms */
} /* namespace elephon */

#include "Algorithms/TrilinearInterpolation.hpp"
#endif /* ELEPHON_ALGORITHMS_TRILINEARINTERPOLATION_H_ */
