/*	This file ReadVASPElephonData.h is part of elephon.
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
 *  Created on: Dec 29, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_IOMETHODS_READVASPELEPHONDATA_H_
#define ELEPHON_IOMETHODS_READVASPELEPHONDATA_H_

#include <string>
#include <vector>

namespace elephon
{
namespace IOMethods
{

class ReadVASPElephonData
{
public:

	void read_data_file(
			std::string const & filename,
			int & regular_dim_x,
			int & regular_dim_y,
			int & regular_dim_z,
			std::vector<double> & regularGridPotential);
};

} /* namespace IOMethods */
} /* namespace elephon */

#endif /* ELEPHON_IOMETHODS_READVASPELEPHONDATA_H_ */
