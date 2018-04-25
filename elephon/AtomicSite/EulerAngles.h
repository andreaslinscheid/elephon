/*	This file EulerAngles.h is part of elephon.
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

#ifndef ELEPHON_ATOMICSITE_EULERANGLES_H_
#define ELEPHON_ATOMICSITE_EULERANGLES_H_

#include <vector>

namespace elephon
{
namespace AtomicSite
{

/** Decompose the 3 Euler angles of a rotation matrix.
 *
 * The convention is z-y-z, as used for the Wigner D matrix
 * See https://www.geometrictools.com/Documentation/EulerAngles.pdf
 * and https://en.wikipedia.org/wiki/Wigner_D-matrix
 *
 * @param[in] rotationMatrix	C ordered unitary rotation matrix, interpreted as the an active rotation matrix
 * @param[out] alpha			angle of the second rotation about the z axis, range -pi <= alpha < pi
 * @param[out] beta				angle of the rotation about the y axis, range 0 <= beta < pi
 * @param[out] gamma			angle of the first rotation about the z axis, range -pi <= gamma < pi
 */
void eulerAngles(
		std::vector<double> const & rotationMatrix,
		double & alpha,
		double & beta,
		double & gamma);

} /* namespace AtomicSite */
} /* namespace elephon */

#endif /* ELEPHON_ATOMICSITE_EULERANGLES_H_ */
