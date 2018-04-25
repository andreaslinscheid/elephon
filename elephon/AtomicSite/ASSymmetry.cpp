/*	This file ASSymmetry.cpp is part of elephon.
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
 *  Created on: Jan 2, 2018
 *      Author: A. Linscheid
 */

#include "AtomicSite/ASSymmetry.h"
#include "AtomicSite/EulerAngles.h"
#include "Algorithms/helperfunctions.hpp"
#include <assert.h>

namespace elephon
{
namespace AtomicSite
{

void
ASSymmetry::initialize(
		int lmax,
		std::vector<double> carthesianSymmetryOperations)
{
	assert(carthesianSymmetryOperations.size() % 9 == 0);
	lMax_ = lmax;
	numSymOps_= carthesianSymmetryOperations.size()/9;
	eulerAngles_.resize(numSymOps_*3);
	rotMatricesPtr_.clear();
	rotMatricesPtr_.reserve(numSymOps_);
	for (int isym = 0 ; isym < numSymOps_; ++isym)
	{
		auto thisMat_ptr = carthesianSymmetryOperations.begin()+isym*9;
		auto rotMat = std::vector<double>(thisMat_ptr, thisMat_ptr + 9);
		bool isProperRotation = Algorithms::helperfunctions::determinant_3by3_matrix(rotMat) > 0.0;
		if (not isProperRotation) // remove the inversion from this matrix, the result is a proper rotation matrix
		{
			for (auto & mij : rotMat )
				mij *= -1.0;
		}
		eulerAngles(rotMat, eulerAngles_[isym*3+0], eulerAngles_[isym*3+1], eulerAngles_[isym*3+2]);

		// active / passive rotation conversion
		auto mapAngle = [] (double & angle) {
			angle = -angle;
			angle -= 2.0*M_PI*std::floor(angle/(2.0*M_PI));
		};
		mapAngle(eulerAngles_[isym*3+1]);

		RadSym rotationOperator;
		rotationOperator.reserve(lMax_+1);
		for (int l = 0; l <= lMax_; ++l)
		{
			WignerDMatrix wd;
			wd.initialize(l, eulerAngles_[isym*3+0], eulerAngles_[isym*3+1], eulerAngles_[isym*3+2], isProperRotation);
			rotationOperator.push_back(std::move(wd));
		}
		rotMatricesPtr_.push_back(std::make_shared<RadSym>(std::move(rotationOperator)));
	}
}

std::shared_ptr<const ASSymmetry::RadSym>
ASSymmetry::get_wigner_rotation_matrices_symop(
		int isymop) const
{
	assert((isymop >= 0)&&(isymop<rotMatricesPtr_.size()));
	return rotMatricesPtr_[isymop];
}

} /* namespace AtomicSite */
} /* namespace elephon */
