/*	This file helperfunctions.hpp is part of elephon.
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
 *  Created on: Oct 7, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_ALGORITHMS_HELPERFUNCTIONS_HPP_
#define ELEPHON_ALGORITHMS_HELPERFUNCTIONS_HPP_

#include "LatticeStructure/Tetrahedron.h"
#include <complex>
#include <vector>
#include <algorithm>

namespace elephon
{
namespace Algorithms
{
namespace helperfunctions
{
template<typename T>
struct cmplxTypeTrait
{
	typedef T realT;
};

template<typename T>
struct cmplxTypeTrait<std::complex<T> >
{
	typedef T realT;
};

template<typename T>
T interpolate_single_cube_realtive(
		double x, double y, double z,
		T f000, T f100, T f110, T f010,
		T f001, T f101, T f111, T f011)
{
	typedef typename cmplxTypeTrait<T>::realT realT;

	//Linear interpolation logic
	auto bilinear_interpol = [] (
			double x, double y,
			T f00, T f10, T f11, T f01)
	{
		T a = f00 * realT(1. - x) + f10 * realT(x);
		T b = f01 * realT(1. - x) + f11 * realT(x);
		return a * realT(1. - y) + b * realT(y);
	};

	auto e = bilinear_interpol(x, y, f000, f100, f110, f010);
	auto f = bilinear_interpol(x, y, f001, f101, f111, f011);
	return e * realT( 1. - z) + f * realT(z);
}

template<typename T>
void gradient_interpolation_single_cube_relative(
			double x, double y, double z,
			T f000, T f100, T f110, T f010,
			T f001, T f101, T f111, T f011,
			T & dx, T & dy, T & dz)
{
	// analytic formula http://paulbourke.net/miscellaneous/interpolation/
	//	V = f000*(1 - x)*(1 - y)*(1 - z) +
	//			f100*x (1 - y)*(1 - z) +
	//			f010*(1 - x)*y*(1 - z) +
	//			f001*(1 - x)*(1 - y)*z +
	//			f101*x*(1 - y)*z +
	//			f011*(1 - x)*y*z +
	//			f110*x*y*(1 - z) +
	//			f111*x*y*z;

	dx= f000*(-1)*(1 - y)*(1 - z) +
		f100*(1 - y)*(1 - z) +
		f010*(-1)*y*(1 - z) +
		f001*(-1)*(1 - y)*z +
		f101*(1 - y)*z +
		f011*(-1)*y*z +
		f110*y*(1 - z) +
		f111*y*z;

	dy= f000*(1 - x)*(-1)*(1 - z) +
		f100*x*(-1)*(1 - z) +
		f010*(1 - x)*(1 - z) +
		f001*(1 - x)*(-1)*z +
		f101*x*(-1)*z +
		f011*(1 - x)*z +
		f110*x*(1 - z) +
		f111*x*z;

	dz= f000*(1 - x)*(1 - y)*(-1) +
		f100*x*(1 - y)*(-1) +
		f010*(1 - x)*y*(-1) +
		f001*(1 - x)*(1 - y) +
		f101*x*(1 - y) +
		f011*(1 - x)*y +
		f110*x*y*(-1) +
		f111*x*y;
}

/** @file Contains simple templated functions that do not fit into a particular category.
 */

template<typename T>
T
triangle_area(std::vector<T> const & p1,
		std::vector<T> const & p2,
		std::vector<T> const & p3)
{
	assert((p1.size() == 3) and (p2.size() == 3) and (p3.size() == 3));
	return 0.5*std::sqrt(std::pow(-(p1[1]*p2[0]) + p1[0]*p2[1] + p1[1]*p3[0] - p2[1]*p3[0] -
		       p1[0]*p3[1] + p2[0]*p3[1],2) +
			std::pow(p1[2]*p2[0] - p1[0]*p2[2] - p1[2]*p3[0] + p2[2]*p3[0] +
		       p1[0]*p3[2] - p2[0]*p3[2],2) +
			   std::pow(-(p1[2]*p2[1]) + p1[1]*p2[2] + p1[2]*p3[1] - p2[2]*p3[1] -
		       p1[1]*p3[2] + p2[1]*p3[2],2));
}

/**
 * Given data values on the 4 corner values, linearly interpolate to points within a tetrahedron.
 *
 * @param ptsInTetra	List of x,y,z coordinates for 1 or more points within a tetrahedron. Layout is x1,y1,z1,x2,...,zN.
 * @param tetra			The tetrahedron object. The coordinates \ref ptsInTetra are transformed into Barycentric coordinates.
 * @param cornerData	Vector of 4 vectors with each an equal size of data values.
 * 						The 4 ith elements are interpolated to the ith element of \ref interpolData
 * @param interpolData
 */
template<typename T>
void
interpolate_within_single_tetrahedron(
		std::vector<double> const & ptsInTetra,
		LatticeStructure::Tetrahedron const & tetra,
		std::vector<std::vector<T>> const & cornerData,
		std::vector<T> & interpolData)
{
	typedef typename cmplxTypeTrait<T>::realT realT;
	assert(ptsInTetra.size()%3 == 0);
	assert( cornerData.size() == 4 );
	int nD = cornerData[0].size();
	assert( (nD == cornerData[1].size()) && (nD == cornerData[2].size()) &&
			(nD == cornerData[3].size()) );

	const int nVThisTetra = ptsInTetra.size()/3;
	std::vector<bool> isInTetra;
	std::vector<double> vectorsBarycentric;
	tetra.check_vectors_inside( ptsInTetra,
								isInTetra,
								vectorsBarycentric);
	assert(std::all_of(isInTetra.begin(), isInTetra.end(), [] (bool a){return a;}));

	interpolData.resize(nD*nVThisTetra);
	for (int ip = 0 ; ip < nVThisTetra ; ++ip)
		for ( int id = 0 ; id < nD ; ++id )
		{
			interpolData[ip*nD+id] =	cornerData[0][id]*realT(vectorsBarycentric[ip*4+0]) +
										cornerData[1][id]*realT(vectorsBarycentric[ip*4+1]) +
										cornerData[2][id]*realT(vectorsBarycentric[ip*4+2]) +
										cornerData[3][id]*realT(vectorsBarycentric[ip*4+3]) ;
			//check for NaN in debug mode
			assert( interpolData[ip*nD+id] == interpolData[ip*nD+id]);
		}
}

inline void
compute_spherical_coords(
		double x, double y, double z,
		double &r, double & theta, double & phi)
{
	r = std::sqrt(std::pow(x,2) + std::pow(y,2) + std::pow(z,2));
	if (r < 1e-10) // by convention, set theta and phi to zero for zero radius
	{
		phi = 0.0;
		theta = 0.0;
		return;
	}
	theta = std::acos(z/r);
	phi = std::atan2(y,x);
}

/**
 * Compute the determinant of a 3x3 matrix.
 *
 * @tparam M			A container that implements the size function and defines the content type as M::value_type
 * @param[in] mat		3x3 matrix in c storage layout
 * @return				the determinant of the matrix
 */
template<class M>
typename M::value_type
determinant_3by3_matrix(M const & mat)
{
	assert(mat.size() == 9);
	typename M::value_type matrix[3][3];
	std::copy(mat.data(), mat.data()+mat.size(), &matrix[0][0]);
	return 	-matrix[0][2]*matrix[1][1]*matrix[2][0] + matrix[0][1]*matrix[1][2]*matrix[2][0] +
			 matrix[0][2]*matrix[1][0]*matrix[2][1] - matrix[0][0]*matrix[1][2]*matrix[2][1] -
			 matrix[0][1]*matrix[1][0]*matrix[2][2] + matrix[0][0]*matrix[1][1]*matrix[2][2];
}

} /* namespace helperfunctions */
} /* namespace Algorithms */
} /* namespace elephon */

#endif /* ELEPHON_ALGORITHMS_HELPERFUNCTIONS_HPP_ */
