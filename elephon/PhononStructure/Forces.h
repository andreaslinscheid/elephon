/*	This file Forces.h is part of elephon.
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
 *  Created on: Feb 23, 2018
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_PHONONSTRUCTURE_FORCES_H_
#define ELEPHON_PHONONSTRUCTURE_FORCES_H_

#include "Auxillary/AlignedVector.h"
#include <memory>

namespace elephon
{
namespace LatticeStructure { class AtomDisplacementCollection; };
namespace IOMethods { class ElectronicStructureCodeInterface; };

namespace PhononStructure
{

class Forces
{
public:

	void initialize(
			std::shared_ptr<const LatticeStructure::AtomDisplacementCollection> displ,
			std::shared_ptr<IOMethods::ElectronicStructureCodeInterface> loader );

	Auxillary::Multi_array<double,2> const & get_forces_for_atom( int atomIndex, int irredDisplIndex) const;

	int get_num_total_irred_displacements() const;

private:

	int totalNumIrredDispl_;

	Auxillary::Multi_array<int,2> indexMap_;

	std::vector<Auxillary::Multi_array<double,2>> forceData_;
};

} /* namespace PhononStructure */
} /* namespace elephon */

#endif /* ELEPHON_PHONONSTRUCTURE_FORCES_H_ */
