/*	This file LatticeModule.cpp is part of elephon.
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
 *  Created on: May 15, 2017
 *      Author: A. Linscheid
 */

#include "LatticeModule.h"
#include <utility>
#include <assert.h>
#include <cmath>

namespace elephon
{
namespace LatticeStructure
{

void LatticeModule::initialize( std::vector<double> latticeMatrix )
{
	latticeMatrix_ = std::move(latticeMatrix);
	assert( latticeMatrix_.size() == 9 );

	//We define alat as the length of the first lattice vector
	//This defines a length scale for Cartesian coordinates
	std::vector<double> a1 = {0.0,0.0,0.0};
	for ( int i = 0 ; i < 3; i++)
		a1[i] = latticeMatrix_[i*3+0];
	alat_ = std::sqrt( a1[0]*a1[0] + a1[1]*a1[1] + a1[2]*a1[2] );
	for (auto &aij : latticeMatrix_)
		aij /= alat_;

	//Compute the reciprocal lattice matrix
	auto a2 = a1;
	auto a3 = a1;
	auto b1 = a1;
	auto b2 = a1;
	auto b3 = a1;
	for ( int i = 0 ; i < 3; i++)
	{
		a1[i] = latticeMatrix_[i*3+0];
		a2[i] = latticeMatrix_[i*3+1];
		a3[i] = latticeMatrix_[i*3+2];
	}
	this->cross_prod(b1,a2,a3);
	this->cross_prod(b2,a3,a1);
	this->cross_prod(b3,a1,a2);

	for ( int i = 0 ; i < 3; i++)
	{
		 reciLatMatrix_[i*3+0] = b1[i];
		 reciLatMatrix_[i*3+1] = b2[i];
		 reciLatMatrix_[i*3+2] = b3[i];
	}
}

std::vector<double> const &
LatticeModule::get_latticeMatrix() const
{
	return latticeMatrix_;
}

std::vector<double> const &
LatticeModule::get_reciprocal_latticeMatrix() const
{
	return reciLatMatrix_;
}

double LatticeModule::get_alat() const
{
	return alat_;
}

void LatticeModule::cross_prod( std::vector<double> const& v1,
		std::vector<double> const& v2,
		std::vector<double> & v1xv2) const
{
	if ( v1xv2.size() != 3 )
		v1xv2 = v1;

	v1xv2[0]=v1[1]*v2[2]-v1[2]*v2[1];
	v1xv2[1]=v1[2]*v2[0]-v1[0]*v2[2];
	v1xv2[2]=v1[0]*v2[1]-v1[1]*v2[0];
}

} /* namespace LatticeStructure */
} /* namespace elephon */
