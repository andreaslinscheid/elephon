/*	This file ElePhonMatModSquare.h is part of elephon.
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
 *  Created on: Nov 16, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_PHONONSTRUCTURE_ELEPHONMATMODSQUARE_H_
#define ELEPHON_PHONONSTRUCTURE_ELEPHONMATMODSQUARE_H_

#include "PhononStructure/ElectronPhononCoupling.h"
#include <memory>
#include <vector>

namespace elephon
{
namespace PhononStructure
{

class ElePhonMatModSquare
{
	void compute_gkkp_mod_square(
			std::shared_ptr<const ElectronPhononCoupling> elphC,
			std::vector<float> & phononModes,
			std::vector<float> & gkkp2) const;
};

} /* namespace PhononStructure */
} /* namespace elephon */

#endif /* ELEPHON_PHONONSTRUCTURE_ELEPHONMATMODSQUARE_H_ */
