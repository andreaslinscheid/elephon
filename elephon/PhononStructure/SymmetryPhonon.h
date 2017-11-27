/*	This file SymmetryPhonon.h is part of elephon.
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
 *  Created on: Nov 26, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_PHONONSTRUCTURE_SYMMETRYPHONON_H_
#define ELEPHON_PHONONSTRUCTURE_SYMMETRYPHONON_H_

#include "LatticeStructure/Symmetry.h"
#include "PhononStructure/Phonon.h"
#include <vector>
#include <memory>

namespace elephon
{
namespace PhononStructure
{

class SymmetryPhonon
{
	void initialize(
			LatticeStructure::Symmetry & sym,
			std::shared_ptr<const PhononStructure::Phonon> ph );

	void calculate_symmetry_q(
			std::vector<double> const & q
			) const;

};

} /* namespace PhononStructure */
} /* namespace elephon */

#endif /* ELEPHON_PHONONSTRUCTURE_SYMMETRYPHONON_H_ */
