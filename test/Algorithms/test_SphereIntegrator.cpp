/*	This file test_SphereIntegrator.cpp is part of elephon.
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
 *  Created on: Jan 31, 2018
 *      Author: A. Linscheid
 */
#include <boost/test/unit_test.hpp>
#include <Algorithms/SphereIntegrator.h>
#include <Auxillary/FunctionTraits.h>
#include <type_traits>
#include <complex>
#include <cmath>

BOOST_AUTO_TEST_SUITE( SphereIntegrator )

struct test_base
{
	double operator() (double a, double b)
	{
		return 1.0;
	}
};

struct test_functor : public test_base {
	mutable bool evaluate_many_has_been_called = false;

	void evaluate_many(double const * theta, double const * phi, int numPts, double * data) const
	{
		evaluate_many_has_been_called = true;
		std::fill(data, data+numPts, 1.0);
	}
};

// randomly pick one integration rule.
inline int choose_rule()
{
	std::vector<int> avail_rules = {3, 5, 7, 9, 11, 13, 15,
			17, 19, 21, 23, 25, 27, 29, 31, 35, 41, 47, 53, 59, 65, 71, 77,
			83, 89, 95, 101, 107, 113, 119, 125, 131 };
	srand(time(nullptr));
	auto ruleIndex = rand() % static_cast<int>(avail_rules.size());
	BOOST_MESSAGE("Selected rule No. " + std::to_string(avail_rules[ruleIndex]));
	return avail_rules[ruleIndex];
}

BOOST_AUTO_TEST_CASE( integrate_sphere_surface )
{
	auto iden = [] (double phi, double theta) {return 1.0;};

	// here come some checks on FunctionTraits - we do this here, because most of the items for
	// testing are here anyway.
	static_assert(std::is_same<elephon::Auxillary::FunctionTraits<decltype(iden)>::result_type, double>::value,
			"Error, return type of lambda not correctly obtained");

	static_assert(std::is_same<elephon::Auxillary::FunctionTraits<test_functor>::result_type, double>::value,
			"Error, return type of inherited operator() not correctly obtained");

	elephon::Algorithms::SphereIntegrator<decltype(iden)> integrator;
	auto unitSphereSurface = integrator.integrate(iden);
	BOOST_CHECK_SMALL((unitSphereSurface- 4.0*M_PI), 1e-8);

	test_functor func;
	elephon::Algorithms::SphereIntegrator<test_functor> integrator2;
	unitSphereSurface = integrator2.integrate(func, choose_rule());
	BOOST_CHECK_SMALL((unitSphereSurface- 4.0*M_PI), 1e-8);
	BOOST_CHECK(func.evaluate_many_has_been_called);
}

struct spherical_harm_Y_l_1_m_0 {
	std::complex<double> operator() (double theta, double phi) const {
		return 0.5*std::sqrt(3.0/M_PI)*std::cos(theta);
	};
};

struct spherical_harm_Y_l_1_m_1 {
	std::complex<double> operator() (double theta, double phi) const {
		return -0.5*std::sqrt(3.0/2.0/M_PI)*std::sin(theta)*
				std::complex<double>(std::cos(phi),std::sin(phi));
	};
};

template<class C1, class C2>
struct prod_functor{
	prod_functor(C1 c1, C2 c2) : c1_(c1), c2_(c2) { };
	std::complex<double> operator() (double theta, double phi) const {
		return c1_(theta,phi)*std::conj(c2_(theta,phi));
	}
	C1 c1_;
	C2 c2_;
};

BOOST_AUTO_TEST_CASE( integrate_spherical_harmonics )
{
	// numerically confirm the normalization of an identical spherical harmonic
	spherical_harm_Y_l_1_m_1 yl1m1;
	prod_functor<spherical_harm_Y_l_1_m_1, spherical_harm_Y_l_1_m_1> f2(yl1m1, yl1m1);
	elephon::Algorithms::SphereIntegrator<decltype(f2)> integrator2;
	auto unitSphereSurface = integrator2.integrate(f2, choose_rule());
	BOOST_CHECK_SMALL(std::abs(unitSphereSurface-std::complex<double>(1.0)), 1e-8); // Ylm's supposed to be normalized

	// numerically confirm the orthogonality of two spherical harmonics
	spherical_harm_Y_l_1_m_0 yl1m0;
	prod_functor<spherical_harm_Y_l_1_m_0, spherical_harm_Y_l_1_m_1> function(yl1m0, yl1m1);

	elephon::Algorithms::SphereIntegrator<decltype(function)> integrator;
	unitSphereSurface = integrator.integrate(function, choose_rule());
	BOOST_CHECK_SMALL(std::abs(unitSphereSurface), 1e-8); // Ylm's supposed to be orthogonal
}

BOOST_AUTO_TEST_SUITE_END()
