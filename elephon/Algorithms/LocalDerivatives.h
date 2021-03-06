/*	This file LocalDerivatives.h is part of elephon.
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
 *  Created on: Sep 28, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_ALGORITHMS_LOCALDERIVATIVES_H_
#define ELEPHON_ALGORITHMS_LOCALDERIVATIVES_H_

#include "LatticeStructure/RegularBareGrid.h"
#include <vector>

namespace elephon
{
namespace Algorithms
{
namespace localDerivatives
{

template<typename T, class DataLoader>
void
compute_derivatives_sqr_polynom(
		int nBnd,
		std::vector<int> const & reducibleKPTIndices,
		std::vector<T> * gradientFieldPtr,
		std::vector<T> * hessianFieldPtr,
		LatticeStructure::RegularBareGrid const & grid,
		DataLoader const & reducibleData);

template<typename T, class DataLoader>
void
compute_derivatives_sqr_polynom_symmetric(
		int nBnd,
		std::vector<int> const & reducibleKPTIndices,
		std::vector<T> * gradientFieldPtr,
		std::vector<T> * hessianFieldPtr,
		LatticeStructure::RegularBareGrid const & grid,
		LatticeStructure::Symmetry const & sym,
		DataLoader const & reducibleData);

} /* namespace localDerivatives */
} /* namespace Algorithms */
} /* namespace elephon */

#include "Algorithms/LocalDerivatives.hpp"
#endif /* ELEPHON_ALGORITHMS_LOCALDERIVATIVES_H_ */
