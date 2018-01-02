/*	This file EulerAngles.cpp is part of elephon.
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

#include "AtomicSite/EulerAngles.h"
#include <assert.h>
#include <cmath>

namespace elephon
{
namespace AtomicSite
{

void eulerAngles(
		std::vector<double> const & rotationMatrix,
		double & alpha,
		double & beta,
		double & gamma)
{
	assert(rotationMatrix.size() == 9);
	auto r02 = rotationMatrix[2];
	auto r10 = rotationMatrix[3]; auto r11 = rotationMatrix[4]; auto r12 = rotationMatrix[5];
	auto r20 = rotationMatrix[6]; auto r21 = rotationMatrix[7]; auto r22 = rotationMatrix[8];

	if (r22 < +1)
	{
		if (r22 > -1)
		{
			alpha = std::atan2(r12, r02);
			beta = std::acos(r22);
			gamma = std::atan2(r21, -r20);
		}
		else // r22 = -1
		{
			// Not a unique solution: gamma - alpha = atan2( r10 , r11 )
			alpha = -std::atan2(r10 , r11);
			beta = M_PI;
			gamma = 0.0;
		}
	}
	else // r22 = +1
	{
		// Not a unique solution: gamma + alpha = atan2( r10 , r11 )
		alpha = std::atan2(r10 , r11);
		beta = 0.0;
		gamma = 0.0;
	}
}

} /* namespace AtomicSite */
} /* namespace elephon */
