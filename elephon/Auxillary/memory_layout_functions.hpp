/*	This file memory_layout_functions.hpp is part of elephon.
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
 *  Created on: Jan 17, 2018
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_AUXILLARY_MEMORY_LAYOUT_FUNCTIONS_HPP_
#define ELEPHON_AUXILLARY_MEMORY_LAYOUT_FUNCTIONS_HPP_

#include <cassert>

namespace elephon
{
namespace Auxillary
{
namespace memlayout
{
/** @file Simple functions defining some layout of multi-dimensional data.
 * @todo Recently found Boost.MultiArray. Look into this ...
 */

/**
 *	Defines the memory layout of angular momentum channels.
 *
 * @param l		The main angular momentum quantum number.
 * @param m		The magnetic quantum number, m is in the range [-l,l]
 * @return		The position of the element (l,m)
 */
inline int
angular_momentum_layout(
		int l,
		int m)
{
	assert((m>=-l)&&(m<=l));
	// sum l'=0,l (2l'+1) = (l-1)l + l = l*l
	return l+m + (l*l);
}

} /* namespace memlayout */
} /* namespace Auxillary */
} /* namespace elephon */

#endif /* ELEPHON_AUXILLARY_MEMORY_LAYOUT_FUNCTIONS_HPP_ */
