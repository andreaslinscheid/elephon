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
#include <boost/multi_array.hpp>
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

/**
 * Wrapper for boost::multi_array to work with linear algebra interface that uses.
 *
 * We redefine several types that have a different meaning in the context of stl containers.
 * In particular size() is not the number of individual elements in the boost multiarray but the number
 * of 'objects' in the first dimension. This makes sense, but it not what we want here.
 * Similarly, value_type is redefined as compared to standard containers.
 */
template<typename T, std::size_t NumDims,
typename Allocator = boost::alignment::aligned_allocator<T, alignedvector::architecture_align>>
class Multi_array : private boost::multi_array<T, NumDims, Allocator>
{
public:
	 using boost::multi_array<T, NumDims, Allocator>::data;
	 using boost::multi_array<T, NumDims, Allocator>::operator[];
	 using boost::multi_array<T, NumDims, Allocator>::begin;
	 using boost::multi_array<T, NumDims, Allocator>::end;
	 using boost::multi_array<T, NumDims, Allocator>::shape;
	 using boost::multi_array<T, NumDims, Allocator>::resize;
	 using boost::multi_array<T, NumDims, Allocator>::reshape;

	 typedef typename boost::multi_array<T, NumDims, Allocator>::element value_type;

	 Multi_array() { };

	 template <class ExtentList>
	 explicit Multi_array(
		      ExtentList const& extents)
	  	  : boost::multi_array<T, NumDims, Allocator>(extents){ };

	typename boost::multi_array<T, NumDims, Allocator>::size_type
	size() const {return boost::multi_array<T, NumDims, Allocator>::num_elements();};
};

} /* namespace Auxillary */
} /* namespace elephon */

#endif /* ELEPHON_AUXILLARY_ALIGNEDVECTOR_H_ */
