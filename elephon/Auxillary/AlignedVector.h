/*	This file AlignedVector.h is part of elephon.
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
 *  Created on: Nov 17, 2017
 *      Author: A. Linscheid
 */

#include <boost/align/aligned_allocator.hpp>
#include <vector>
#include <complex>

#ifndef ELEPHON_AUXILLARY_ALIGNEDVECTOR_H_
#define ELEPHON_AUXILLARY_ALIGNEDVECTOR_H_

namespace elephon
{
namespace Auxillary
{
namespace alignedvector
{
const int architecture_align = 128;

template<class T, std::size_t Alignment = architecture_align>
using aligned_vector = std::vector<T,
    boost::alignment::aligned_allocator<T, Alignment> >;

typedef aligned_vector<char, architecture_align> BV;

typedef aligned_vector<float, architecture_align> FV;

typedef aligned_vector<double, architecture_align> DV;

typedef aligned_vector<std::complex<float>, architecture_align> CV;

typedef aligned_vector<std::complex<double>, architecture_align> ZV;

inline bool
is_aligned(const void * ptr, std::uintptr_t alignment) noexcept
{
    auto iptr = reinterpret_cast<std::uintptr_t>(ptr);
    return !(iptr % alignment);
}

template<typename T>
bool check_data_is_aligned(T * ptr)
{
	return is_aligned(reinterpret_cast<void*>(ptr), architecture_align);
}

} /* namespace alignedvector */
} /* namespace Auxillary */
} /* namespace elephon */

#endif /* ELEPHON_AUXILLARY_ALIGNEDVECTOR_H_ */
