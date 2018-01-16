/*	This file SymmetryOperation.cpp is part of elephon.
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
#include <cmath>

namespace elephon
{
namespace symmetry
{

SymmetryOperation::SymmetryOperation(
		std::vector<int>::const_iterator pointGroupLatticeBasisBegin,
		std::vector<int>::const_iterator pointGroupLatticeBasisEnd,
		std::vector<double>::const_iterator fractLatticeBasisBegin,
		std::vector<double>::const_iterator fractLatticeBasisEnd,
		std::vector<double>::const_iterator pointGroupCarthBasisBegin,
		std::vector<double>::const_iterator pointGroupCarthBasisEnd,
		std::vector<double>::const_iterator fractCarthBasisBegin,
		std::vector<double>::const_iterator fractCarthBasisEnd )
{
	assert(std::distance(pointGroupLatticeBasisBegin, pointGroupLatticeBasisEnd) == 9);
	assert(std::distance(fractLatticeBasisBegin, fractLatticeBasisEnd) == 3);
	assert(std::distance(pointGroupCarthBasisBegin, pointGroupCarthBasisEnd) == 9);
	assert(std::distance(fractCarthBasisBegin, fractCarthBasisEnd) == 3);

	std::copy(pointGroupLatticeBasisBegin, pointGroupLatticeBasisEnd, ptgroup);
	std::copy(fractLatticeBasisBegin, fractLatticeBasisEnd, fracTrans);
	std::copy(pointGroupCarthBasisBegin, pointGroupCarthBasisEnd, ptgCart);
	std::copy(fractCarthBasisBegin, fractCarthBasisEnd, fracTransCart);
}

void
SymmetryOperation::rotate( std::vector<int> & v ) const
{
	assert(v.size()==3);
	auto b = v;
	for ( int i = 0; i < 3; ++i)
		v[i] = ptgroup[i*3+0]*b[0]+ptgroup[i*3+1]*b[1]+ptgroup[i*3+2]*b[2];
}

void
SymmetryOperation::apply( std::vector<double> & v, bool latticePeriodic ) const
{
	assert(v.size()==3);
	auto b = v;
	for ( int i = 0; i < 3; ++i)
	{
		v[i] = ptgroup[i*3+0]*b[0]+ptgroup[i*3+1]*b[1]+ptgroup[i*3+2]*b[2]+fracTrans[i];
		if ( latticePeriodic )
			v[i] -= std::floor(v[i]+0.5);
	}
};

void
SymmetryOperation::rotate_cart( std::vector<double> & v ) const
{
	assert(v.size()==3);
	auto b = v;
	for ( int i = 0; i < 3; ++i)
		v[i] = ptgCart[i*3+0]*b[0]+ptgCart[i*3+1]*b[1]+ptgCart[i*3+2]*b[2];
};


} /* namespace symmetry */
} /* namespace elephon */
