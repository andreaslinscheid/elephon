/*	This file Eliashberg_helperfunction.hpp is part of elephon.
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
 *  Created on: May 18, 2018
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_ELIASHBERGEQUATIONS_ELIASHBERG_HELPERFUNCTION_HPP_
#define ELEPHON_ELIASHBERGEQUATIONS_ELIASHBERG_HELPERFUNCTION_HPP_

#include "EliashbergEquations/EliashbergModule.h"
#include <cassert>
#include <cmath>

namespace elephon {
namespace EliashbergEquations {

inline double
fermi_matsubara_frequency_of_index(int index, int numMats, double inverseTemperature)
{
	assert(numMats%2 == 0);
	return M_PI / inverseTemperature * ( 2*(index - numMats/2) + 1);
}

inline double
bose_matsubara_frequency_of_index(int index, int numMats, double inverseTemperature)
{
	assert(numMats%2 == 1);
	return 2.0*M_PI / inverseTemperature * (index - numMats/2);
}

inline int
fermi_number_matsubara_frequencies(double temperature, double MatsubaraCutoff)
{
	double prefactor = M_PI / Algorithms::helperfunctions::inverse_temperature_eV(temperature);;
	return 2*std::floor((MatsubaraCutoff-prefactor)/(2.0*prefactor));
}

} /* namespace EliashbergEquations */
} /* namespace elephon */

#endif /* ELEPHON_ELIASHBERGEQUATIONS_ELIASHBERG_HELPERFUNCTION_HPP_ */
