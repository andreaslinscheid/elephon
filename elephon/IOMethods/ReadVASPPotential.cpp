/*	This file ReadVASPPotential.cpp is part of elephon.
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
 *  Created on: Apr 26, 2018
 *      Author: A. Linscheid
 */

#include "IOMethods/ReadVASPPotential.h"

namespace elephon
{
namespace IOMethods
{

void
ReadVASPPotential::set_filepath(std::string filepath)
{
	filepath_ = std::move(filepath);
}

std::string const &
ReadVASPPotential::get_default_filename() const
{
	return defaultFileName_;
}

} /* namespace IOMethods */
} /* namespace elephon */
