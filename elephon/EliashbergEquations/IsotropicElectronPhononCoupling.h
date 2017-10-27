/*	This file IsotropicElectronPhononCoupling.h is part of elephon.
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
 *  Created on: Oct 25, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_ELIASHBERGEQUATIONS_ISOTROPICELECTRONPHONONCOUPLING_H_
#define ELEPHON_ELIASHBERGEQUATIONS_ISOTROPICELECTRONPHONONCOUPLING_H_

#include "EliashbergEquations/MatsubaraBaseBoson.h"
#include <string>
#include <vector>

namespace elephon
{
namespace EliashbergEquations
{

class IsotropicElectronPhononCoupling : public MatsubaraBaseBoson
{
public:

	void initialize(
			std::string const & filename,
			double temperature,
			double energyCutoff);

	void initialize(
			std::vector<double> const & frequencies,
			std::vector<double> const & a2F,
			double temperature,
			double energyCutoff);

};

} /* namespace EliashbergEquations */
} /* namespace elephon */

#endif /* ELEPHON_ELIASHBERGEQUATIONS_ISOTROPICELECTRONPHONONCOUPLING_H_ */
