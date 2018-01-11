/*	This file test_CubeSplineInterpolation.cpp is part of elephon.
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
 *  Created on: Jan 11, 2018
 *      Author: A. Linscheid
 */
#include <boost/test/unit_test.hpp>
#include "Algorithms/CubeSplineInterpolation.h"
#include "Auxillary/CubicPolynomial.h"
#include "Algorithms/helperfunctions.hpp"
#include <vector>
#include <complex>

BOOST_AUTO_TEST_SUITE( CubeSplineInterpolation )

template<typename TMesh, typename TData>
void
test_setup()
{
	typedef typename elephon::Algorithms::helperfunctions::cmplxTypeTrait<TData>::realT RT;

	//Initialize with spline test data that can be solved analytically
	//	(see http://en.wikipedia.org/wiki/Spline_interpolation)
	//This test data leads to two polynomials that are supposed to have the	derivatives
	//	(k0 = -0.6875 , k1 = -0.1250) and (k0 = -0.1250 , k1 = 1.5625 )
	std::vector<TMesh> xValues;
	std::vector<TData> dataSet;
	xValues.push_back(-1.0);xValues.push_back(0.0);xValues.push_back(3.0);
	dataSet.push_back( 0.5);dataSet.push_back(0.0);dataSet.push_back(3.0);

	elephon::Algorithms::CubeSplineInterpolation<TMesh, TData> cubeSpline;
	cubeSpline.initialize(xValues, dataSet);

	elephon::Auxillary::CubicPolynomial<TMesh, TData> firstPolynom(xValues[0],xValues[1],
			dataSet[0],dataSet[1],-0.6875,-0.1250);
	elephon::Auxillary::CubicPolynomial<TMesh, TData> secondPolynom(xValues[1],xValues[2],
			dataSet[1],dataSet[2],-0.1250,1.5625);

	const int numEvals = 100;
	for ( int i = 0 ; i < numEvals ; i++)
	{
		// sample the interval
		TMesh x = cubeSpline.min_range() + ((cubeSpline.max_range() - cubeSpline.min_range() )/numEvals)*i;

		TData value;
		RT diffValue;
		if ( x < xValues[1])
		{
			firstPolynom.evaluate(x, value);
			diffValue = std::abs(value - cubeSpline(x));
		}
		else
		{
			secondPolynom.evaluate(x, value);
			diffValue = std::abs(value - cubeSpline(x));
		}
		BOOST_CHECK_SMALL(diffValue, RT(5e-5) );
	}
};

BOOST_AUTO_TEST_CASE( CubeSplineInterpolation_analytic_example )
{
	// check all allowed combinations for the CubeSplineInterpolation
	test_setup<double, double>();
	test_setup<double, std::complex<double>>();
	test_setup<float, double>();
	test_setup<float, std::complex<double>>();
	test_setup<double, float>();
	test_setup<double, std::complex<float>>();
	test_setup<float, float>();
	test_setup<float, std::complex<float>>();
}

BOOST_AUTO_TEST_SUITE_END()
