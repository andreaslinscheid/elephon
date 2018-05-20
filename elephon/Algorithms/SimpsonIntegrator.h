/*	This file SimpsonIntegrator.h is part of elephon.
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
 *  Created on: May 17, 2018
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_ALGORITHMS_SIMPSONINTEGRATOR_H_
#define ELEPHON_ALGORITHMS_SIMPSONINTEGRATOR_H_

#include "Auxillary/FunctionTraits.h"
#include <iterator>

namespace elephon {
namespace Algorithms {

template<class D>
struct SimpsonIntegrator
{
	typedef typename std::iterator_traits<typename D::const_iterator>::value_type T;

	T integrate(D const & dataArray, T stepWidth);

	T integrate(D const & dataArray, D const & positions);

	T integrate(typename D::const_iterator dataArrayBegin, typename D::const_iterator dataArrayEnd,
			typename D::const_iterator positionsBegin, typename D::const_iterator positionsEnd);
};

} /* namespace Algorithms */
} /* namespace elephon */

#include "Algorithms/SimpsonIntegrator.hpp"
#endif /* ELEPHON_ALGORITHMS_SIMPSONINTEGRATOR_H_ */
