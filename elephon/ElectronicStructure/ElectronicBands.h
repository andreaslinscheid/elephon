/*	This file ElectronicBands.h is part of elephon.
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
 *  Created on: May 20, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_ELECTRONICSTRUCTURE_ELECTRONICBANDS_H_
#define ELEPHON_ELECTRONICSTRUCTURE_ELECTRONICBANDS_H_

#include "LatticeStructure/RegularGrid.h"
#include <vector>

namespace elephon
{
namespace ElectronicStructure
{

class ElectronicBands
{
public:

	void initialize(
			std::vector<double> const & kpoints,
			int numBands,
			std::vector<double> bandData,
			LatticeStructure::RegularGrid grid);

	int get_nBnd() const;

	void generate_reducible_grid_bands(
			std::vector<int> const & bIndices,
			std::vector<double> & bands) const;

	LatticeStructure::RegularGrid const & get_grid() const;
private:

	int nBnd_;

	std::vector<double> dataIrred_;

	LatticeStructure::RegularGrid grid_;
};

} /* namespace ElectronicStructure */
} /* namespace elephon */

#endif /* ELEPHON_ELECTRONICSTRUCTURE_ELECTRONICBANDS_H_ */
