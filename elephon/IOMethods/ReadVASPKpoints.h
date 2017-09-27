/*	This file ReadVASPKpoints.h is part of elephon.
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
 *  Created on: Sep 27, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_IOMETHODS_READVASPKPOINTS_H_
#define ELEPHON_IOMETHODS_READVASPKPOINTS_H_

#include "LatticeStructure/LatticeModule.h"
#include <string>

namespace elephon
{
namespace IOMethods
{

class ReadVASPKpoints
{
public:

	void read_kpoints(
			std::string filename,
			LatticeStructure::LatticeModule const & lattice);

	std::vector<double> const & get_grid_shift() const;

	std::vector<int> const & get_grid_dim() const;

private:

	std::string filename_;

	std::vector<int> dim_ = {0, 0, 0};

	std::vector<double> gridShift_ = {0.0, 0.0, 0.0};
};

} /* namespace IOMethods */
} /* namespace elephon */

#endif /* ELEPHON_IOMETHODS_READVASPKPOINTS_H_ */
