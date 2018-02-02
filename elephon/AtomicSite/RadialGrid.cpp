/*	This file RadialGrid.cpp is part of elephon.
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

#include "AtomicSite/RadialGrid.h"
#include "Algorithms/CubeSplineInterpolation.h"

namespace elephon
{
namespace AtomicSite
{

void
RadialGrid::initialize(
		std::vector<double> center,
		double radius,
		std::vector<double> points)
{
	numR_ = points.size();
	center_ = std::move(center);
	radius_ = radius;
	radialPoints_ = std::move(points);
	assert(center_.size() == 3);
}

int
RadialGrid::get_num_R() const
{
	return numR_;
}

double
RadialGrid::get_radius(int iR) const
{
	assert((iR >= 0)&&(iR<numR_));
	return radialPoints_[iR];
}

std::vector<double> const &
RadialGrid::get_center() const
{
	assert(center_.size() == 3);
	return center_;
}

void
RadialGrid::transform(symmetry::SymmetryOperation const & sop)
{
	sop.apply(center_);
}

} /* namespace AtomicSite */
} /* namespace elephon */
