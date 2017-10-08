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

#include <complex>

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


} /* namespace helperfunctions */
} /* namespace Algorithms */
} /* namespace elephon */

#endif /* ELEPHON_ALGORITHMS_HELPERFUNCTIONS_HPP_ */
