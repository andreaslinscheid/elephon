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
#include "Algorithms/helperfunctions.hpp"
#include <utility>
#include <assert.h>
#include <cmath>
#include <stdexcept>
#include <iostream>

namespace elephon
{
namespace LatticeStructure
{

LatticeModule::LatticeModule()
{
	this->initialize( {1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0} );
}

LatticeModule::LatticeModule(std::vector<double> latticeMatrix )
{
	this->initialize( std::move(latticeMatrix) );
}

void
LatticeModule::initialize( std::vector<double> latticeMatrix )
{
	latticeMatrix_ = std::move(latticeMatrix);
	assert( latticeMatrix_.size() == 9 );

	std::vector<double> a1 = this->get_lattice_vector(0);
	auto a2 = this->get_lattice_vector(1);
	auto a3 = this->get_lattice_vector(2);

	//Compute the reciprocal lattice matrix
	auto b1 = a1;
	auto b2 = a1;
	auto b3 = a1;

	Algorithms::helperfunctions::cross_prod(a2,a3,b1);
	Algorithms::helperfunctions::cross_prod(a3,a1,b2);
	Algorithms::helperfunctions::cross_prod(a1,a2,b3);

	vol_ = a1[0]*b1[0]+a1[1]*b1[1]+a1[2]*b1[2];
	assert( vol_ > 0 );

	reciLatMatrix_.resize(9);
	for ( int i = 0 ; i < 3; i++)
	{
		 reciLatMatrix_[i*3+0] = b1[i]/vol_;
		 reciLatMatrix_[i*3+1] = b2[i]/vol_;
		 reciLatMatrix_[i*3+2] = b3[i]/vol_;
	}

	//We define alat as the length of the first lattice vector
	//This defines a length scale for Cartesian coordinates
	//We measure the reciprocal basis in units of 2pi/alat
	alat_ = std::sqrt( a1[0]*a1[0] + a1[1]*a1[1] + a1[2]*a1[2] );
	for (auto &aij : latticeMatrix_)
		aij /= alat_;
	for (auto &bij : reciLatMatrix_)
		bij *= alat_;

	std::vector<double> b(9,0.0);
	for ( int i = 0 ; i < 3 ; ++i)
		for ( int j = 0 ; j < 3 ; ++j)
			for ( int k = 0 ; k < 3 ; ++k)
				b[i*3+j] += reciLatMatrix_[k*3+i]*latticeMatrix_[k*3+j];
	std::vector<double> c(9,0.0);
	for ( int i = 0 ; i < 3 ; ++i)
		for ( int j = 0 ; j < 3 ; ++j)
			for ( int k = 0 ; k < 3 ; ++k)
				c[i*3+j] += latticeMatrix_[k*3+i]*reciLatMatrix_[k*3+j];
	for ( int i = 0 ; i < 3 ; ++i)
		for ( int j = 0 ; j < 3 ; ++j)
			if ( (std::abs(b[i*3+j] - (i==j?1.0:0.0)) > 1e-6) or (std::abs(c[i*3+j] - (i==j?1.0:0.0)) > 1e-6) )
				throw std::runtime_error("Problem of non-orthogonal lattice and reciprocal matrix");
}

LatticeModule
LatticeModule::build_supercell(std::vector<int> const & supercellDim) const
{
	if(supercellDim.size() == 3)
	{
		LatticeStructure::LatticeModule newLattice;
		std::vector<double> newLatticeMatrix(9, 0.0);
		for ( int i = 0 ; i < 3; ++i)
			for ( int j = 0 ; j < 3; ++j)
				newLatticeMatrix[i*3+j] = this->get_latticeMatrix()[i*3+j]*supercellDim[i]*this->get_alat();
		newLattice.initialize(newLatticeMatrix);
		return newLattice;
	}
	Auxillary::Multi_array<int,2> matrix(boost::extents[3][3]);
	assert(supercellDim.size() == 9);
	std::copy(supercellDim.begin(),supercellDim.end(), matrix.data());
	return this->build_supercell(matrix);
}

LatticeModule
LatticeModule::build_supercell(Auxillary::Multi_array<int,2> const & supercellMatrix) const
{
	LatticeStructure::LatticeModule newLattice;
	std::vector<double> newLatticeMatrix(9, 0.0);
	assert((supercellMatrix.shape()[0] == 3) && (supercellMatrix.shape()[1] == 3));
	// arrange loop such that lattice vector i gets multiplied from the left with the supercell matrix.
	for ( int i = 0 ; i < 3; ++i)
		for ( int j = 0 ; j < 3; ++j)
			for ( int k = 0 ; k < 3; ++k)
				newLatticeMatrix[i*3+j] += supercellMatrix[j][k]*this->get_latticeMatrix()[i*3+k]*this->get_alat();
	newLattice.initialize(newLatticeMatrix);
	return newLattice;
}

//
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

double
LatticeModule::get_volume() const
{
	return vol_;
}

double
LatticeModule::get_reci_volume() const
{
	return std::pow(2*M_PI,3)/vol_;
}

std::vector<double>
LatticeModule::get_lattice_vector(int n) const
{
	assert( (n >= 0) && (n < 3) && (latticeMatrix_.size() == 9));
	return {latticeMatrix_[0*3+n],latticeMatrix_[1*3+n],latticeMatrix_[2*3+n]};
}

std::vector<double>
LatticeModule::get_reci_lattice_vector(int n) const
{
	assert( (n >= 0) && (n < 3) && (reciLatMatrix_.size() == 9));
	return {reciLatMatrix_[0*3+n],reciLatMatrix_[1*3+n],reciLatMatrix_[2*3+n]};
}

void
LatticeModule::direct_to_cartesian_matrix(double * mat, int nelem) const
{
	assert( (nelem%9 == 0) && ( latticeMatrix_.size() == 9 ) && (reciLatMatrix_.size() == 9));
	std::vector<double> b(9);
	int nc = nelem/9;
	for ( int ic = 0 ; ic < nc; ic++)
	{
		std::fill(b.data(),b.data()+9, 0.0 );
		for ( int i = 0 ; i < 3 ; ++i)
			for ( int j = 0 ; j < 3 ; ++j)
				for ( int k = 0 ; k < 3 ; ++k)
					for ( int l = 0 ; l < 3 ; ++l)
						b[i*3+j] += latticeMatrix_[i*3+k]*mat[ic*9+k*3+l]*reciLatMatrix_[j*3+l];
		std::copy(b.data(),b.data()+9,&mat[ic*9]);
	}
}

void
LatticeModule::reci_cartesian_to_direct_matrix(double * mat, int nelem) const
{
	assert( (nelem%9 == 0) && ( latticeMatrix_.size() == 9 ) && (reciLatMatrix_.size() == 9));
	std::vector<double> b(9);
	int nc = nelem/9;
	for ( int ic = 0 ; ic < nc; ic++)
	{
		std::fill(b.data(),b.data()+9, 0.0 );
		for ( int i = 0 ; i < 3 ; ++i)
			for ( int j = 0 ; j < 3 ; ++j)
				for ( int k = 0 ; k < 3 ; ++k)
					for ( int l = 0 ; l < 3 ; ++l)
						b[i*3+j] += latticeMatrix_[k*3+i]*mat[ic*9+k*3+l]*reciLatMatrix_[l*3+j];
		std::copy(b.data(),b.data()+9,&mat[ic*9]);
	}
}

void
LatticeModule::reci_direct_to_cartesian(std::vector<double> & v) const
{
	this->reci_direct_to_cartesian(v.data(),v.size());
}

void
LatticeModule::reci_direct_to_cartesian(double * v, int nelem) const
{
	assert( (nelem%3 == 0) && ( reciLatMatrix_.size() == 9 ));
	double b[3];
	int nc = nelem/3;
	for ( int ic = 0 ; ic < nc; ic++)
	{
		std::copy(&v[ic*3],&v[ic*3]+3,b);
		for ( int i = 0 ; i < 3; i++)
			v[ic*3+i] = reciLatMatrix_[i*3+0]*b[0]+reciLatMatrix_[i*3+1]*b[1]+reciLatMatrix_[i*3+2]*b[2];
	}
}

void
LatticeModule::reci_direct_to_cartesian_2pibya(std::vector<double> & v) const
{
	this->reci_direct_to_cartesian(v.data(),v.size());
	for ( auto &vi : v )
		vi *= 2.0*M_PI/alat_;
}

void
LatticeModule::reci_direct_to_cartesian_2pibya(double * p , int nelem) const
{
	this->reci_direct_to_cartesian(p,nelem);
	for ( int i = 0 ; i < nelem; ++i )
		p[i] *= 2.0*M_PI/alat_;
}


void
LatticeModule::direct_to_cartesian(double p[3]) const
{
	this->direct_to_cartesian(p,3);
}

void
LatticeModule::direct_to_cartesian(double * p, int nelem) const
{
	assert( latticeMatrix_.size() == 9 );
	assert(nelem%3 == 0);
	double b[3];
	for (int iele = 0 ; iele < nelem/3; ++iele)
	{
		double * p_this_ele = p + 3*iele;
		std::copy(p_this_ele,p_this_ele+3,b);
		for ( int i = 0 ; i < 3; i++)
			p_this_ele[i] = latticeMatrix_[i*3+0]*b[0]+latticeMatrix_[i*3+1]*b[1]+latticeMatrix_[i*3+2]*b[2];
	}
}

void LatticeModule::direct_to_cartesian(std::vector<double> & v) const
{
	assert( (v.size()%3 == 0) && ( latticeMatrix_.size() == 9 ));
	double b[3];
	int nc = int(v.size())/3;
	for ( int ic = 0 ; ic < nc; ic++)
	{
		std::copy(&v[ic*3],&v[ic*3]+3,b);
		for ( int i = 0 ; i < 3; i++)
			v[ic*3+i] = latticeMatrix_[i*3+0]*b[0]+latticeMatrix_[i*3+1]*b[1]+latticeMatrix_[i*3+2]*b[2];
	}
}

void
LatticeModule::direct_to_cartesian_angstroem(double * p, int nelem) const
{
	this->direct_to_cartesian(p,nelem);
	for (int i = 0 ; i < nelem; ++i)
		p[i] *= alat_;
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
	assert( (v.size()%3 == 0) && ( reciLatMatrix_.size() == 9 ));
	double b[3];
	int nc = int(v.size())/3;
	for ( int ic = 0 ; ic < nc; ic++)
	{
		std::copy(&v[ic*3],&v[ic*3]+3,b);
		for ( int i = 0 ; i < 3; i++)
			v[ic*3+i] = reciLatMatrix_[0*3+i]*b[0]+reciLatMatrix_[1*3+i]*b[1]+reciLatMatrix_[2*3+i]*b[2];
	}
}

} /* namespace LatticeStructure */
} /* namespace elephon */
