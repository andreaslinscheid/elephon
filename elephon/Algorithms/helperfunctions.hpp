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

#include "Auxillary/AlignedVector.h"
#include "Auxillary/UnitConversion.h"
#include <complex>
#include <vector>
#include <algorithm>
#include <utility>
#include <type_traits>
#include <assert.h>
#include <typeinfo>
#include <fstream>

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
template<typename T, class Tetrahedron>
void
interpolate_within_single_tetrahedron(
		std::vector<double> const & ptsInTetra,
		Tetrahedron const & tetra,
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

namespace detail {
// a type trait that calls operator[] if it exists or operator() else
template<class T>
struct has_op_bracket
{
	template<class A, class B>
	static auto check(A * p, B * q) -> typename std::is_same<decltype((*p)(std::declval<int>())), decltype((*q)(std::declval<int>()))>::type;

	template<typename  ... Args>
	static auto check(Args ... args) -> std::false_type;

	static constexpr bool value = decltype(check(static_cast<T *>(0),static_cast<T *>(0)))::value;
};

template<class T, bool has_op_bracket> struct map_trait { };
template<class T>
struct map_trait<T,true> {
	template<typename ...Args>
	static int call(T const & a, Args ... args) {return a(args...);};
};
template<class T>
struct map_trait<T,false> {
	static int call(T const & a, int index) {return a[index];};
};
}


/**
 * Apply the inverse of the rotation operator to the scalar field data at a given point in 3D space.
 *
 * Since the symmetry operation rotates a position vector, after a transformation, the data at
 * a given point belongs to the target of the transformation. Hence for a forward transformation
 * the data field transform with the inverse.
 *
 * @param[in,out] scalarFieldBegin	Iterator pointing to the beginning of the range.
 * @param[in,out] scalarFieldEnd	Iterator pointing to the element behind the last of the range.
 * @param[in] map					A class that implements the operator[] and tells for a given index i where i ends up
 */
template<class iterator, class indexMap>
void transform_scalar_field_cart(
		iterator scalarFieldBegin, iterator scalarFieldEnd,
		indexMap const & map)
{
	assert((std::distance(scalarFieldBegin, scalarFieldEnd)>0));
	const int nD = std::distance(scalarFieldBegin, scalarFieldEnd);
	typedef typename std::iterator_traits<iterator>::value_type T;
	auto runnerIt = scalarFieldBegin;

	Auxillary::alignedvector::aligned_vector<T> bufferField(scalarFieldBegin, scalarFieldEnd);
	// re-shuffle the scalar field
	for (int iD = 0 ; iD < nD; ++iD)
	{
		int mapped_index = detail::map_trait<indexMap,detail::has_op_bracket<indexMap>::value>::call(map,iD);
		assert((mapped_index>=0) && (mapped_index < bufferField.size()));
		bufferField[mapped_index] = *runnerIt++;
	}
	std::copy(bufferField.begin(), bufferField.end(), scalarFieldBegin);
}

/**
 * Apply the inverse of the rotation operator to the vector field data at a given point in 3D space.
 *
 * Since the symmetry operation rotates a position vector, after a transformation, the data at
 * a given point belongs to the target of the transformation. Hence for a forward transformation
 * the data field transform with the inverse.
 *
 * @param[in,out] vectorFieldBegin	Iterator pointing to the beginning of the range. 3 consqutive elements
 * 									are field value in x,y,z respectively. std::distance(vectorFieldBegin, vectorFieldEnd) >= 3 required.
 * 									The range must be evenly devidible by 3.
 * @param[in,out] vectorFieldEnd	Iterator pointing to the element behind the last of the range.
 * 									std::distance(vectorFieldBegin, vectorFieldEnd) >= 3 required.
 * 									The range must be evenly devidible by 3.
 * @param[in] mat					A type with consecutive storage that holds the rotation matrix in c-layout.
 * @param[in] map					A class that implements the operator[] and tells for a given index i where i ends up
 */
template<class iterator, class M, class indexMap>
void transform_vector_field_cart(
		iterator vectorFieldBegin, iterator vectorFieldEnd,
		M const & mat,
		indexMap const & map)
{
	assert(mat.size() == 9);
	typename M::value_type matrix[3][3];
	std::copy(mat.data(), mat.data()+mat.size(), &matrix[0][0]);
	assert((std::distance(vectorFieldBegin, vectorFieldEnd)%3 == 0)&&(std::distance(vectorFieldBegin, vectorFieldEnd)>0));
	const int nD = std::distance(vectorFieldBegin, vectorFieldEnd)/3;
	typedef typename std::iterator_traits<iterator>::value_type T;
	auto currentIterator = vectorFieldBegin;

	T LocalBuffer[3];
	Auxillary::alignedvector::aligned_vector<T> bufferField(vectorFieldBegin, vectorFieldEnd);
	// re-shuffle the vector field and rotate the vector at each point
	for (int iD = 0 ; iD < nD; ++iD)
	{
		LocalBuffer[0] = *currentIterator++;
		LocalBuffer[1] = *currentIterator++;
		LocalBuffer[2] = *currentIterator++;
		int mapped_index = detail::map_trait<indexMap,detail::has_op_bracket<indexMap>::value>::call(map,iD);
		assert((mapped_index>=0) && (mapped_index < bufferField.size()/3));
		bufferField[mapped_index*3+0] = T(matrix[0][0])*LocalBuffer[0] + T(matrix[0][1])*LocalBuffer[1] + T(matrix[0][2])*LocalBuffer[2];
		bufferField[mapped_index*3+1] = T(matrix[1][0])*LocalBuffer[0] + T(matrix[1][1])*LocalBuffer[1] + T(matrix[1][2])*LocalBuffer[2];
		bufferField[mapped_index*3+2] = T(matrix[2][0])*LocalBuffer[0] + T(matrix[2][1])*LocalBuffer[1] + T(matrix[2][2])*LocalBuffer[2];
	}

	std::copy(bufferField.begin(), bufferField.end(), vectorFieldBegin);
}

/**
 * Compute the integer closest to a floating point number.
 *
 * @tparam a floating point type
 *
 * @param[in] v	the floating point number
 * @return	the integer closest to the floating point number.
 */
template<typename T>
int nint(T v)
{
	return std::floor(v+T(0.5));
}

template<class VT, typename ReadT = typename VT::value_type>
void
read_binary_file(const char * filename_ctr, VT & data)
{
	std::ifstream file(filename_ctr, std::ios::binary);
	if ( ! file.good() )
		throw std::runtime_error(std::string("Failed to open binary file ") + filename_ctr + " for reading!");

	file.seekg(0, std::ios::beg);
	std::int64_t begin = file.tellg();
	file.seekg(0, std::ios::end);
	std::int64_t end = file.tellg();
	std::int64_t size = end - begin;

	if (size % sizeof(ReadT) != 0)
		throw std::runtime_error(std::string("Error reading ") + filename_ctr
				+ ": interpreting content as " + typeid(ReadT).name() + " leads to left-over bytes");

	const std::int64_t numElem = size/sizeof(ReadT);
	data.resize(numElem);
	file.seekg(0, std::ios::beg);

	if ( typeid(ReadT) == typeid(typename VT::value_type))
	{
		file.read(reinterpret_cast<char*>(data.data()), size);
	}
	else
	{
		ReadT  * buffer = new ReadT [numElem];
		file.read(reinterpret_cast<char*>(buffer), size);
		std::copy(buffer, buffer+numElem, data.data());
		delete [] buffer;
	}
}

/**
 * Compute the cross product of two vectors.
 *
 * @param[in] v1		First vector. Must be of size 3.
 * @param[in] v2		Second vector. Must be of size 3.
 * @param[out] v1xv2	Result vector, resized to 3 if its size is not 3.
 */
template<class V>
void cross_prod( V const& v1, V const& v2, V & v1xv2)
{
	assert(v1.size() == 3);
	assert(v2.size() == 3);
	if ( v1xv2.size() != 3 )
		v1xv2 = v1;

	v1xv2[0]=v1[1]*v2[2]-v1[2]*v2[1];
	v1xv2[1]=v1[2]*v2[0]-v1[0]*v2[2];
	v1xv2[2]=v1[0]*v2[1]-v1[1]*v2[0];
}

template<typename T>
T inverse_temperature_eV(T temperature)
{
	assert(temperature > T(0));
	return T(1.0)/Auxillary::units::BOLTZMANN_CONSTANT_IN_EV_PER_K/temperature;
}

} /* namespace helperfunctions */
} /* namespace Algorithms */
} /* namespace elephon */

#endif /* ELEPHON_ALGORITHMS_HELPERFUNCTIONS_HPP_ */
