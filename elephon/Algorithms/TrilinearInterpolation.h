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
#include "LatticeStructure/RegularBareGrid.h"

namespace elephon
{
namespace Algorithms
{

class TrilinearInterpolation
{
public:
	TrilinearInterpolation(
			LatticeStructure::RegularBareGrid grid);

	/**
	 * requiredGridIndices will be given in a consecutively ordered grid index
	 * where where x is the slowest running: (ix,iy,iz) <==> (ix*grid[y]+iy)*grid[z]+iz
	 */
	void data_query(
			std::vector<double> listOfPoints,
			std::vector<int> & requiredGridIndices);
	/**
	 *
	 */
	void interpolate(
			int nDataPtsPerPoint,
			std::vector<double> const& gridDataForRequiredIndices,
			std::vector<double> & pointsData) const;

	template<typename T>
	void interpolate_within_single_cube(
			std::vector<double> const & ptsInCube,
			std::vector<std::vector<T>> const & cornerData,
			std::vector<T> & interpolData) const;

private:

	LatticeStructure::RegularBareGrid  grid_;

	std::vector<double> listOfPoints_;

	std::vector<LatticeStructure::RegularBareGrid::GridCube> usedGridCubes_;

	std::vector<int> conseqPtsRegularGrid_;

	void get_cell_vectors( std::vector<double> const& k,
			std::vector<double> & lowerCorner,
			std::vector<double> & upperCorner ) const;

	template<typename T>
	void interpolate(
			int nDataPtsPerPoint,
			std::vector<T> const& gridDataForRequiredIndices,
			std::vector<T> & pointsData) const;

};

} /* namespace Algorithms */
} /* namespace elephon */

#include "Algorithms/TrilinearInterpolation.hpp"
#endif /* ELEPHON_ALGORITHMS_TRILINEARINTERPOLATION_H_ */
