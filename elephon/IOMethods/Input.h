/*	This file Input.h is part of elephon.
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
 *  Created on: Apr 24, 2017
 *      Author: A. Linscheid
 */

#ifndef IOMETHODS_INPUT_H_
#define IOMETHODS_INPUT_H_

#include <boost/program_options.hpp>

namespace elephon {
namespace IOMethods {

class Input {
public:
	Input( int argc, char const * const* argv );

	std::vector<size_t> get_super_cell_size() const;

	size_t get_numFS() const;
private:

	boost::program_options::variables_map vm_;
};

} /* namespace IOMethods */
} /* namespace elephon */

#endif /* IOMETHODS_INPUT_H_ */
