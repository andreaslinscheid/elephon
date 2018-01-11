/*	This file test_CubicPolynomial.cpp is part of elephon.
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
#include "Auxillary/CubicPolynomial.h"
#include "Algorithms/helperfunctions.hpp"
#include <vector>
#include <complex>

BOOST_AUTO_TEST_SUITE( CubicPolynomial )

// the main test logic is templated, since we want to check for (complex) float, double independently
template<typename TMesh, typename TData>
void
test_setup()
{
	typedef typename elephon::Algorithms::helperfunctions::cmplxTypeTrait<TData>::realT RT;

	// test for a constant polynomial that noting weird happens
	const TData valueOfPolynom = TData(2.0);
	elephon::Auxillary::CubicPolynomial<TMesh, TData> cpd_constant(
			TMesh(1.0), TMesh(2.0),	// range
			valueOfPolynom, valueOfPolynom,	// value
			TData(0.0), TData(0.0));			// derivatives
	TData value;
	cpd_constant.evaluate(1.5, value);
	BOOST_CHECK_SMALL( std::fabs(value - valueOfPolynom), RT(5e-5) );

	// some real analytic test data
	std::vector<double> xValuesD = {-0.989218429384234234,-0.2319384234234,
			-0.0012134234876,0.62347289433245,1.3453158523413245};
	std::vector<double> dataArrayD = {-5.3891789678248282332,.8209848419213039104,
			0.9997192497522201019,1.9744090986534953532,12.0114397015586715940};
	std::vector<double> derivativeDataD = {17.8055531218330275340,1.5949988948483670073,
			0.2327484514620032387,4.9287797622237035407,25.6058553748127175706};

	std::vector<TMesh> xValues(xValuesD.begin(), xValuesD.end());
	std::vector<TData> dataArray(dataArrayD.begin(), dataArrayD.end());
	std::vector<TData> derivativeData(derivativeDataD.begin(), derivativeDataD.end());

	// check that the function evaluates to the input reference values
	elephon::Auxillary::CubicPolynomial<TMesh, TData> cpd_testData(
			xValues.front(), xValues.back()+std::numeric_limits<TMesh>::epsilon(),
			dataArray.front(), dataArray.back(),
			derivativeData.front(), derivativeData.back());

	for ( int ixp = 0 ; ixp < xValues.size() ; ixp++)
	{
		TData value, derivative;
		cpd_testData.evaluate(xValues[ixp], value);
		BOOST_CHECK_SMALL( std::fabs(value - dataArray[ixp]), RT(5e-5) );
		cpd_testData.evaluate_derivative(xValues[ixp], derivative);
		BOOST_CHECK_SMALL( std::fabs(derivative - derivativeData[ixp]), RT(5e-5) );
	}
};

BOOST_AUTO_TEST_CASE( CubicPolynomial_analytic_example )
{
	// check all allowed combinations for the polynomial
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
