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

template<class iterator>
void
Symmetry::rotate_cartesian(
		int isym,
		iterator fieldCartBegin,
		iterator fieldCartEnd) const
{
	assert( std::distance(fieldCartBegin,fieldCartEnd)%3 == 0 );
	assert( isym < numRotations_ );

	typename std::iterator_traits<iterator>::value_type buff[3];

	int numComponents = std::distance(fieldCartBegin,fieldCartEnd)/3;
	for ( int ic = 0; ic < numComponents; ++ic)
	{
		std::copy(fieldCartBegin+ic*3,fieldCartBegin+(ic+1)*3, &buff[0]);
		for ( int xi = 0; xi < 3; ++xi)
		{
			*(fieldCartBegin+ic*3+xi) = symmetriesCartesian_[(isym*3+xi)*3+0]*buff[0]
							+symmetriesCartesian_[(isym*3+xi)*3+1]*buff[1]
							+symmetriesCartesian_[(isym*3+xi)*3+2]*buff[2];
		}
	}
}

template<class iterator>
void
Symmetry::rotate_matrix_cartesian(int isym,
		iterator matrixFieldCartBegin,
		iterator matrixFieldCartEnd) const
{
	int elem = std::distance(matrixFieldCartBegin,matrixFieldCartEnd);
	assert( elem%9 == 0 );
	assert( isym < numRotations_ );

	std::vector<double> b(9);

	int numComponents = elem/9;
	for ( int ic = 0; ic < numComponents; ++ic)
	{
		std::copy(matrixFieldCartBegin+ic*9,matrixFieldCartBegin+(ic+1)*9,std::begin(b));
		std::fill(matrixFieldCartBegin+ic*9,matrixFieldCartBegin+(ic+1)*9, 0.0 );
		auto G = &symmetriesCartesian_[isym*9];
		for ( int i = 0; i < 3; ++i)
			for ( int j = 0; j < 3; ++j)
				for ( int k = 0; k < 3; ++k)
					for ( int l = 0; l < 3; ++l)
						*(matrixFieldCartBegin+ic*9+i*3+l) += G[i*3+j]*b[j*3+k]*G[l*3+k];
	}
}


} /* namespace LatticeStructure */
} /* namespace elephon */
