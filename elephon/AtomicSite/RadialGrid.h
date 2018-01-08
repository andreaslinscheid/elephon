/*	This file RadialGrid.h is part of elephon.
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
 *  Created on: Jan 4, 2018
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_ATOMICSITE_RADIALGRID_H_
#define ELEPHON_ATOMICSITE_RADIALGRID_H_

#include <vector>

namespace elephon
{
namespace AtomicSite
{

class RadialGrid
{
public:

	void initialize(
			std::vector<double> center,
			double radius,
			std::vector<double> points);

	int get_num_R() const;

	/**
	 * Perform a cubic spline interpolation of data on the grid.
	 *
	 * @param rValues				List of points where the data should be interpolated to.
	 * @param nDataPerRValue		Number of blocks of data of size {rValues}, each of which will be
	 * 								interpolated individually.
	 * @param gridDataBegin			Constant iterator to the beginning of the grid data.
	 * 								 The data to be interpolated is [gridDataBegin, gridDataEnd[
	 * @param gridDataEnd			Constant iterator to end of the grid data.
	 * 								 The data to be interpolated is [gridDataBegin, gridDataEnd[
	 * @param interpolDataBegin		Iterator to the range where the interpolated data should be placed.
	 */
	template<class ConstIterator, class Iterator>
	void interpolate(
			std::vector<double> const & rValues,
			int nDataPerRValue,
			ConstIterator gridDataBegin,
			ConstIterator gridDataEnd,
			Iterator interpolDataBegin) const;
private:

	int numR_ = 1;

	std::vector<double> center_ = {0.0, 0.0, 0.0};

	double radius_ = 0;

	std::vector<double> radialPoints_ = {0.0};
};

} /* namespace AtomicSite */
} /* namespace elephon */

#include "AtomicSite/RadialGrid.hpp"
#endif /* ELEPHON_ATOMICSITE_RADIALGRID_H_ */
