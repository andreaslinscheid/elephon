/*	This file PhononGrid.h is part of elephon.
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
 *  Created on: Oct 2, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_PHONONSTRUCTURE_PHONONGRID_H_
#define ELEPHON_PHONONSTRUCTURE_PHONONGRID_H_

#include "PhononStructure/Phonon.h"
#include "LatticeStructure/DataRegularGrid.h"
#include "LatticeStructure/RegularSymmetricGrid.h"
#include "LatticeStructure/TetrahedraGrid.h"
#include <vector>
#include <string>

namespace elephon
{
namespace PhononStructure
{

class PhononGrid : public LatticeStructure::DataRegularGrid<double>
{
public:

	int num_modes() const;

	void write_phonon_dos_file(std::string const & filename,
			std::vector<double> frequencies = std::vector<double>(),
			std::shared_ptr<const LatticeStructure::TetrahedraGrid> tetra = nullptr,
			std::shared_ptr<const Phonon> ph = nullptr) const;

	void write_phonon_tetrahedra_dos_file(
			std::string const & filename,
			std::vector<double> frequencies) const;

private:
};

} /* namespace PhononStructure */
} /* namespace elephon */

#endif /* ELEPHON_PHONONSTRUCTURE_PHONONGRID_H_ */
