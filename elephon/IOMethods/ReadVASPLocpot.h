/*	This file ReadVASPPotcar.h is part of elephon.
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
 *  Created on: May 14, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_IOMETHODS_READVASPLOCPOT_H_
#define ELEPHON_IOMETHODS_READVASPLOCPOT_H_

#include <string>
#include <vector>

namespace elephon
{
namespace IOMethods
{

class ReadVASPLocpot
{
public:

	void read_scf_potential(std::string const & filename,
			std::vector<int> & dims,
			std::vector<double> & potential) const;
private:

};

} /* namespace IOMethods */
} /* namespace elephon */

#endif /* ELEPHON_IOMETHODS_READVASPLOCPOT_H_ */
