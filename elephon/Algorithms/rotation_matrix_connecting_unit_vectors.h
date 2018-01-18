/*	This file rotation_matrix_connecting_unit_vectors.h is part of elephon.
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
 *  Created on: Jan 18, 2018
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_ALGORITHMS_ROTATION_MATRIX_CONNECTING_UNIT_VECTORS_H_
#define ELEPHON_ALGORITHMS_ROTATION_MATRIX_CONNECTING_UNIT_VECTORS_H_

#include <boost/multi_array.hpp>
#include <vector>

namespace elephon
{
namespace Algorithms
{
/** @file */

/**
 *	Given two unit vectors 1 and 2 this method computes a rotation matrix that takes 1 into 2.
 *
 *	The convention is a multiplication from the left, such that R*r1 = r2.
 *
 * @param[in] unitVectorOrig	3 element vector with norm 1 that is taken to unitVectorResult with matrix multiplication from the left.
 * @param[in] unitVectorResult	3 element vector with norm 1 that is the result of the unitVectorOrig multiplied by the result.
 * @param[out] rotationMatrix	assigned the matrix R that takes unitVectorOrig to unitVectorResult by a left multiplication.
 */
inline void
rotation_matrix_connecting_unit_vectors(
		std::vector<double> const & unitVectorOrig,
		std::vector<double> const & unitVectorResult,
		boost::multi_array<double,2> & rotationMatrix);

} /* namespace Algorithms */
} /* namespace elephon */

#endif /* ELEPHON_ALGORITHMS_ROTATION_MATRIX_CONNECTING_UNIT_VECTORS_H_ */
