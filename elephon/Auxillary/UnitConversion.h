/*	This file UnitConversion.h is part of elephon.
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
 *  Created on: Nov 1, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_AUXILLARY_UNITCONVERSION_H_
#define ELEPHON_AUXILLARY_UNITCONVERSION_H_

namespace elephon
{
namespace Auxillary
{
namespace units
{

/**
 * Convert the matrix inverse of the band structure hession to an effective mass
 * in units of the bare electron rest mass.
 *
 * Internal units of the code are eV for energy and Angstroems for length leading to
 * a band structure hessian in units of [ d^2 E / dki dkj ] =  eV / A^2.
 * The mass tensor is M^-1 = 1 / hbar^2 [ d^2 E / dki dkj ] leading to the conversion
 * below to obtain M from [ d^2 E / dki dkj ]^-1 in units of me, the electron rest mass.
 */
const double INVERSE_EV_TIMES_A2_TO_ME = 7.61996389393771;

} /* namespace units */
} /* namespace Auxillary */
} /* namespace elephon */

#endif /* ELEPHON_AUXILLARY_UNITCONVERSION_H_ */
