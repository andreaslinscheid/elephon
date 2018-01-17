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
#include "Auxillary/memory_layout_functions.hpp"
#include "Algorithms/LinearAlgebraInterface.h"
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

template<class VT>
void
SymmetryOperation::rotate_radial_data(
		int lMax,
		int numDataPerLM,
		VT & dataToBeTransformed) const
{
	Algorithms::LinearAlgebraInterface linalg;
	assert(radSym_ptr_->size() >= lMax);

	VT buffer(numDataPerLM*(2*lMax+1));
	for (int l = 0 ; l <= lMax ; ++l)
	{
		int nl = 2*l+1;
		auto const & wignerD = (*radSym_ptr_)[l];
		assert(wignerD.view_as_matrix().size() == nl*nl);
		auto wignerD_ptr = wignerD.view_as_matrix().data();
		auto expansionData_ptr = dataToBeTransformed.data() +
				Auxillary::memlayout::angular_momentum_layout(l,-l)*numDataPerLM;
		linalg.call_gemm('t', 'n',	// See declaration comments; the Wigner matrix is transposed because we transform the coefficients, while the
									//		operator acts on the functions.
				nl, numDataPerLM, nl, // W is a (2*l+1) x (2*l+1) matrix, the expansion data we interpret as a (2*l+1) x numDataPerLM matrix.
				decltype(*wignerD_ptr)(1.0), wignerD_ptr, nl,
				expansionData_ptr, numDataPerLM,
				decltype(*wignerD_ptr)(0.0), buffer.data(), numDataPerLM );

		std::copy(buffer.begin(), buffer.begin()+numDataPerLM*nl, expansionData_ptr);
	}
}

} /* namespace symmetry */
} /* namespace elephon */
