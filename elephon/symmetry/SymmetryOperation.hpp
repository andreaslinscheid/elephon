/*	This file SymmetryOperation.hpp is part of elephon.
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
 *  Created on: Jan 16, 2018
 *      Author: A. Linscheid
 */

#include "symmetry/SymmetryOperation.h"
#include <cassert>

namespace elephon
{
namespace symmetry
{

template<class VT>
void
SymmetryOperation::rotate_matrix_cart(VT & m) const
{
	typedef typename VT::value_type T;
	assert(m.size()==9);
	auto buffer = m;
	std::fill(buffer.begin(), buffer.end(), T(0));
	for ( int i = 0; i < 3; ++i)
		for ( int j = 0; j < 3; ++j)
			for ( int k = 0; k < 3; ++k)
				for ( int l = 0; l < 3; ++l)
					buffer[i*3+j] += ptgCart[i*3+k]*m[k*3+l]*ptgCart[j*3+l];
}

inline double
SymmetryOperation::get_carth_rot_matrix(int i, int j) const
{
	assert((i>=0)&&(i<3));
	assert((j>=0)&&(j<3));
	return ptgCart[i*3+j];
}

inline double
SymmetryOperation::get_carth_frac_trans(int i) const
{
	assert((i>=0)&&(i<3));
	return fracTransCart[i];
}

inline int
SymmetryOperation::get_lat_rot_matrix (int i, int j) const
{
	assert((i>=0)&&(i<3));
	assert((j>=0)&&(j<3));
	return ptgroup[i*3+j];
}

inline int
SymmetryOperation::get_lat_frac_trans(int i) const
{
	assert((i>=0)&&(i<3));
	return fracTrans[i];
}

} /* namespace symmetry */
} /* namespace elephon */
