/*	This file rotation_matrix_connecting_unit_vectors.cpp is part of elephon.
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
 *  Created on: Jan 18, 2018
 *      Author: A. Linscheid
 */

#include "Algorithms/rotation_matrix_connecting_unit_vectors.h"
#include <cassert>

namespace elephon
{
namespace Algorithms
{

void
rotation_matrix_connecting_unit_vectors(
		std::vector<double> const & unitVectorOrig,
		std::vector<double> const & unitVectorResult,
		boost::multi_array<double,2> & rotationMatrix)
{
	// formula according to https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
	assert(unitVectorOrig.size() == 3);
	assert(unitVectorResult.size() == 3);
	rotationMatrix.resize(boost::extents[3][3]);

	auto cross_prod = [] (std::vector<double> const & a, std::vector<double> const & b) {
		return std::vector<double>{ 	a[1]*b[2] - a[2]*b[1],
										a[2]*b[0] - a[0]*b[2],
										a[0]*b[1] - a[1]*b[0]};
	};
	auto dprod = [] (std::vector<double> const & a, std::vector<double> const & b){
		return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
	};
	assert(std::abs(dprod(unitVectorOrig,unitVectorOrig) - 1.0) < 1e-6);
	assert(std::abs(dprod(unitVectorResult,unitVectorResult) - 1.0) < 1e-6);

	// set to id
	std::fill_n(rotationMatrix.data(), rotationMatrix.num_elements(), 0.0);
	for (int i = 0 ; i < 3; ++i)
		rotationMatrix[i][i] = 1.0;

	std::vector<double> const & a = unitVectorOrig;
	std::vector<double> const & b = unitVectorResult;

	std::vector<double> v = cross_prod(a, b);
	double s = std::sqrt(dprod(v,v));
	double c = dprod(a,b);
	// a || b
	if (std::abs(s) < 1e-8)
	{
		if ( c > 0 ) // a equals b
			return; // identity matrix

		if ( c < 0 ) // a equals -b, pick any axis orthogonal to a and rotate by pi
		{
			// any vector perpendicular to a satisfies v*a == 0, v=(x,y,z)
			// for simplicity lets take x = 0, or y =0 and require that
			// the vector v is norm 1
			double x,y,z;
			if ( std::pow(a[1],2) < std::pow(a[0],2) )
			{
				// a is more parallel to x than y, in this case we take y == 0
				x = a[2]/std::sqrt(std::pow(a[0],2) + std::pow(a[2],2));
				y = 0.0;
				z = -a[0]/std::sqrt(std::pow(a[0],2) + std::pow(a[2],2));
			}
			else
			{
				// take x = 0;
				x = 0;
				y = a[2]/std::sqrt(std::pow(a[1],2) + std::pow(a[2],2));
				z = -a[1]/std::sqrt(std::pow(a[1],2) + std::pow(a[2],2));
			}

			double norm = std::sqrt(x*x+y*y+z*z);
			x /= norm;
			y /= norm;
			z /= norm;

			// define a rotation matrix of pi about the axis (x,y,z), assuming |x,y,z| == 1
			std::vector<double> rotationPi{std::pow(x,2) - std::pow(y,2) - std::pow(z,2), 2*x*y, 2*x*z,
					   2*x*y, -std::pow(x,2) + std::pow(y,2) - std::pow(z,2), 2*y*z,
					   2*x*z, 2*y*z, -std::pow(x,2) - std::pow(y,2) + std::pow(z,2)};

			// this is also our rotation matrix, copy it to the variable
			std::copy(rotationPi.begin(), rotationPi.end(), rotationMatrix.data());
			return;
		}
	}

	// R = I + vx + vx^2 (1-c)/s
	std::vector<double> skew_matrix_v{0, -v[2], v[1], v[2], 0, -v[0], -v[1], v[0], 0};
	std::vector<double> skew_matrix_v_square{-std::pow(v[1],2) - std::pow(v[2],2),v[0]*v[1],v[0]*v[2],
			   v[0]*v[1],-std::pow(v[0],2) - std::pow(v[2],2),v[1]*v[2],
			   v[0]*v[2],v[1]*v[2],-std::pow(v[0],2) - std::pow(v[1],2)};

	for (int i = 0 ; i < 3 ; ++i)
		for (int j = 0 ; j < 3 ; ++j)
			rotationMatrix[i][j] += skew_matrix_v[i*3+j] + (1.0 - c)/(s*s)*skew_matrix_v_square[i*3+j];
}

} /* namespace Algorithms */
} /* namespace elephon */
