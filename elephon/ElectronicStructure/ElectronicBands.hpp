/*	This file ElectronicBands.hpp is part of elephon.
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
 *  Created on: Sep 25, 2017
 *      Author: A. Linscheid
 */

#include "ElectronicStructure/ElectronicBands.h"
#include "Algorithms/LocalDerivatives.h"

namespace elephon
{
namespace ElectronicStructure
{

template<typename T>
void
ElectronicBands::compute_derivatives_sqr_polynom(
		std::vector<int> const & bandIndices,
		std::vector<int> const & reducibleKPTIndices,
		std::vector<T> * gradientFieldPtr,
		std::vector<T> * hessianFieldPtr ) const
{
	// this lambda take the job of mapping a regular reducible k grid index and a "local" band index
	// into the irreducible zone and the full band-context band index.
	auto translate_bands = [&] (int ikr, int ib) {
		return (*this)( grid_.get_maps_red_to_irreducible()[ikr], bandIndices[ib]);
	};

	Algorithms::localDerivatives::compute_derivatives_sqr_polynom<T>(
			bandIndices.size(),
			reducibleKPTIndices,
			gradientFieldPtr,
			hessianFieldPtr,
			grid_.view_bare_grid(),
			translate_bands );
}

} /* namespace ElectronicStructure */
} /* namespace elephon */
