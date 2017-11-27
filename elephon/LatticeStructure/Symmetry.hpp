/*	This file Symmetry.hpp is part of elephon.
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
 *  Created on: Jul 1, 2017
 *      Author: A. Linscheid
 */

#include "LatticeStructure/Symmetry.h"
#include <algorithm>

namespace elephon
{
namespace LatticeStructure
{

template<typename T>
void
Symmetry::rotate(
		int isym,
		typename std::vector<T>::iterator fieldDirectSpaceBegin,
		typename std::vector<T>::iterator fieldDirectSpaceEnd,
		bool latticePeriodic) 								const
{
	assert( std::distance(fieldDirectSpaceBegin,fieldDirectSpaceEnd)%3 == 0 );
	assert( isym < this->get_num_symmetries() );
	auto fb = fieldDirectSpaceBegin; //shorter names for clarity
	auto fe = fieldDirectSpaceEnd;

	T buff[3];

	int numComponents = std::distance(fb,fe)/3;
	for ( int ic = 0; ic < numComponents; ++ic)
	{
		auto fc = fb+ic*3;
		std::copy(fc,fc+3,buff);
		if ( not isReciprocalSpace_ )
			for ( int xi = 0; xi < 3; ++xi)
			{
				*(fc+xi) = symmetries_[(isym*3+xi)*3+0]*buff[0]
					      +symmetries_[(isym*3+xi)*3+1]*buff[1]
						  +symmetries_[(isym*3+xi)*3+2]*buff[2];
				if ( latticePeriodic )
					*(fc+xi) -= std::floor( *(fc+xi)+0.5);
			}
		else
			for ( int xi = 0; xi < 3; ++xi)
			{
				*(fc+xi) = symmetriesReciprocal_[(isym*3+xi)*3+0]*buff[0]
						  +symmetriesReciprocal_[(isym*3+xi)*3+1]*buff[1]
						  +symmetriesReciprocal_[(isym*3+xi)*3+2]*buff[2];
				if ( latticePeriodic )
					*(fc+xi) -= std::floor( *(fc+xi)+0.5);
			}
	}
}

template<class VT>
void
Symmetry::Sop::rotate_matrix_cart(VT & m) const
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

} /* namespace LatticeStructure */
} /* namespace elephon */
