/*	This file ReadVASPSymmetries.h is part of elephon.
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
 *  Created on: May 16, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_IOMETHODS_READVASPSYMMETRIES_H_
#define ELEPHON_IOMETHODS_READVASPSYMMETRIES_H_

#include <string>
#include <vector>

namespace elephon
{
namespace IOMethods
{

class ReadVASPSymmetries
{
public:

	void read_file(std::string filename );

	std::vector<int> const& get_symmetries() const;

	std::vector<double> const& get_fractionTranslations() const;

	bool get_time_revesal_symmetry() const;

private:

	bool timeReversal_ = false;

	std::vector<int> symmetries_;

	std::vector<double> fractionTranslations_;

	void parse_symmetry_blocks(std::string const & fcontent,
			std::vector<std::string> & blocks) const;
};

} /* namespace IOMethods */
} /* namespace elephon */

#endif /* ELEPHON_IOMETHODS_READVASPSYMMETRIES_H_ */
