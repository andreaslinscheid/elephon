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

	//The basis vectors
	std::vector<double> a1 = {0.0,0.0,0.0};
	auto a2 = a1;
	auto a3 = a1;
	auto b1 = a1;
	auto b2 = a1;
	auto b3 = a1;

	//Compute the reciprocal lattice matrix
	for ( int i = 0 ; i < 3; i++)
	{
		a1[i] = latticeMatrix_[i*3+0];
		a2[i] = latticeMatrix_[i*3+1];
		a3[i] = latticeMatrix_[i*3+2];
	}
	this->cross_prod(a2,a3,b1);
	this->cross_prod(a3,a1,b2);
	this->cross_prod(a1,a2,b3);

	double unitCellVolume = a1[0]*b1[0]+a1[1]*b1[1]+a1[2]*b1[2];
	assert( unitCellVolume > 0 );

	reciLatMatrix_.resize(9);
	for ( int i = 0 ; i < 3; i++)
	{
		 reciLatMatrix_[0*3+i] = b1[i]/unitCellVolume;
		 reciLatMatrix_[1*3+i] = b2[i]/unitCellVolume;
		 reciLatMatrix_[2*3+i] = b3[i]/unitCellVolume;
	}

	//We define alat as the length of the first lattice vector
	//This defines a length scale for Cartesian coordinates
	//We measure the reciprocal basis in units of 2pi/alat
	alat_ = std::sqrt( a1[0]*a1[0] + a1[1]*a1[1] + a1[2]*a1[2] );
	for (auto &aij : latticeMatrix_)
		aij /= alat_;
	for (auto &bij : reciLatMatrix_)
		bij *= alat_;

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

void LatticeModule::direct_to_cartesian(std::vector<double> & v) const
{
	assert( (v.size()%3 == 0) && ( latticeMatrix_.size() == 9 ));
	std::vector<double> b(3);
	int nc = int(v.size())/3;
	for ( int ic = 0 ; ic < nc; ic++)
	{
		std::copy(&v[ic*3],&v[ic*3]+3,b.data());
		for ( int i = 0 ; i < 3; i++)
			v[ic*3+i] = latticeMatrix_[i*3+0]*b[0]+latticeMatrix_[i*3+1]*b[1]+latticeMatrix_[i*3+2]*b[2];
	}
}

void
LatticeModule::direct_to_cartesian_angstroem(std::vector<double> & v) const
{
	this->direct_to_cartesian(v);
	for ( auto &xi : v )
		xi *= alat_;
}

void LatticeModule::cartesian_to_direct(std::vector<double> & v) const
{
	assert( (v.size() == 3) && ( reciLatMatrix_.size() == 9 ));
	auto b = v;
	for ( int i = 0 ; i < 3; i++)
		v[i] = reciLatMatrix_[i*3+0]*b[0]+reciLatMatrix_[i*3+1]*b[1]+reciLatMatrix_[i*3+2]*b[2];
}

} /* namespace LatticeStructure */
} /* namespace elephon */
