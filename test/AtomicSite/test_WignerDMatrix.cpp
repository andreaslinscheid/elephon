/*	This file test_WignerDMatrix.cpp is part of elephon.
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
#include "AtomicSite/WignerDMatrix.h"
#include <cstdlib>
#include <ctime>

BOOST_AUTO_TEST_SUITE( WignerDMatrix )

BOOST_AUTO_TEST_CASE( test_explicit_example )
{
	// This lambda tests the examples for l=0, l=1 and l=2 against
	// analytically known forms of the D matrix
	auto test_wigner_D_angles = [] (double alpha, double beta, double gamma)
	{

		elephon::AtomicSite::WignerDMatrix D0;
		D0.initialize(0, alpha, beta, gamma, /* proper rotatation = */true);
		BOOST_CHECK_CLOSE(std::real(D0(0,0)), 1.0, 1e-10); // scalar
		BOOST_CHECK_CLOSE(std::imag(D0(0,0)), 0.0, 1e-10);

		auto matrix_layout = [] (int l, int m, int mp) {
			return (l+m)*(2*l+1)+(l+mp);
		};

		auto phaseFactor = [&] (int m, int mp) {
			return std::complex<double>(std::cos(alpha*m), std::sin(alpha*m))
					*std::complex<double>(std::cos(gamma*mp), std::sin(gamma*mp));
		};

		int l = 1;
		elephon::AtomicSite::WignerDMatrix D1;
		D1.initialize(l, alpha, beta, gamma, /* proper rotatation = */ true);
		std::vector<std::complex<double>> result(3*3);

		// examples from Wikipedia [https://en.wikipedia.org/wiki/Wigner_D-matrix]
		result[matrix_layout(l,-1,-1)] 	= phaseFactor(-1,-1)*((1.0+std::cos(beta))/2.0            );
		result[matrix_layout(l,-1, 0)] 	= phaseFactor(-1, 0)*( -std::sin(beta)/std::sqrt(2.0)     );
		result[matrix_layout(l,-1, 1)]	= phaseFactor(-1, 1)*((1.0-std::cos(beta))/2.0            );
		result[matrix_layout(l, 0,-1)] 	= phaseFactor( 0,-1)*( std::sin(beta)/std::sqrt(2.0)      );
		result[matrix_layout(l, 0, 0)] 	= phaseFactor( 0, 0)*(std::cos(beta)                      );
		result[matrix_layout(l, 0, 1)] 	= phaseFactor( 0, 1)*( -std::sin(beta)/std::sqrt(2.0)     );
		result[matrix_layout(l, 1,-1)] 	= phaseFactor( 1,-1)*((1.0-std::cos(beta))/2.0            );
		result[matrix_layout(l, 1, 0)] 	= phaseFactor( 1, 0)*( std::sin(beta)/std::sqrt(2.0)      );
		result[matrix_layout(l, 1, 1)] 	= phaseFactor( 1, 1)*((1.0+std::cos(beta))/2.0            );
		for (int m = -l ; m <= l; ++m )
			for (int mp = -l ; mp <= l; ++mp )
			{
				BOOST_CHECK_SMALL(std::abs(D1(m,mp)-result[matrix_layout(l,m,mp)]), 1e-10); // vector
			}

		// examples from mathematica. I love computer algebra ...
		std::vector<double> smallD = { std::pow(std::cos(beta/2.),4),
				 -2*std::pow(std::cos(beta/2.),3)*std::sin(beta/2.),
				 (std::sqrt(1.5)*std::pow(std::sin(beta),2))/2.,
				 -(std::pow(std::sin(beta/2.),2)*std::sin(beta)),
				 std::pow(std::sin(beta/2.),4),
				 ((1 + std::cos(beta))*std::sin(beta))/2.,
				 (std::cos(beta) + std::cos(2*beta))/2.,
				 -(std::sqrt(1.5)*std::cos(beta)*std::sin(beta)),
				   (std::cos(beta) - std::cos(2*beta))/2.,
				   -(std::pow(std::sin(beta/2.),2)*std::sin(beta)),
				   (std::sqrt(1.5)*std::pow(std::sin(beta),2))/2.,
				   std::sqrt(1.5)*std::cos(beta)*std::sin(beta),(1 + 3*std::cos(2*beta))/4.,
				   -(std::sqrt(1.5)*std::cos(beta)*std::sin(beta)),
				   (std::sqrt(1.5)*std::pow(std::sin(beta),2))/2.,
				   std::pow(std::sin(beta/2.),2)*std::sin(beta),
				   (std::cos(beta) - std::cos(2*beta))/2.,
				   std::sqrt(1.5)*std::cos(beta)*std::sin(beta),
				   (std::cos(beta) + std::cos(2*beta))/2.,
				   -2*std::pow(std::cos(beta/2.),3)*std::sin(beta/2.),std::pow(std::sin(beta/2.),4),
				   std::pow(std::sin(beta/2.),2)*std::sin(beta),(std::sqrt(1.5)*std::pow(std::sin(beta),2))/2.,
				   ((1 + std::cos(beta))*std::sin(beta))/2.,std::pow(std::cos(beta/2.),4)};

		l = 2;
		assert(smallD.size() == (2*l+1)*(2*l+1));
		elephon::AtomicSite::WignerDMatrix D2;
		D2.initialize(l, alpha, beta, gamma, /* proper rotatation = */true);
		for (int m = -l ; m <= l; ++m )
			for (int mp = -l ; mp <= l; ++mp )
			{
				BOOST_CHECK_SMALL(std::abs(D2(m,mp)-phaseFactor(m, mp)*smallD[matrix_layout(l,m,mp)]), 1e-10); // rank 2 tensor
			}
	}; // test_wigner_D_angles

	// test all branches for the angle beta
	// the angles alpha and gamma are complex phases and less prone to errors
	// ( (beta > 0) && (beta <= M_PI/2.0) )
	test_wigner_D_angles(/*alpha = */0.0,
						 /*beta = */ M_PI/2,
						 /*gamma = */0.0 );

	// ( beta == 0 )
	test_wigner_D_angles(/*alpha = */0.0,
						 /*beta = */ 0.0,
						 /*gamma = */0.0 );

	// ( (beta > M_PI/2.0) && (beta < M_PI))
	test_wigner_D_angles(/*alpha = */0.0,
						 /*beta = */ M_PI/2.0+0.1,
						 /*gamma = */0.0 );

	// ( beta == M_PI)
	test_wigner_D_angles(/*alpha = */0.0,
						 /*beta = */ M_PI,
						 /*gamma = */0.0 );

	// ( (beta > M_PI) && (beta < 3.0*M_PI/2.0))
	test_wigner_D_angles(/*alpha = */0.0,
						 /*beta = */ M_PI+0.1,
						 /*gamma = */0.0 );

	// ( (beta >= 3.0*M_PI/2.0) && (beta < 2.0*M_PI))
	test_wigner_D_angles(/*alpha = */0.0,
						 /*beta = */ 3.0*M_PI/2.0+0.1,
						 /*gamma = */0.0 );

	std::srand(std::time(NULL));

	// generate 3 random angles in the range [0,2pi]
	double alpha = 2.0*M_PI*static_cast<double>(std::rand())/static_cast<double>(RAND_MAX);
	double beta = 2.0*M_PI*static_cast<double>(std::rand())/static_cast<double>(RAND_MAX);
	double gamma = 2.0*M_PI*static_cast<double>(std::rand())/static_cast<double>(RAND_MAX);

	BOOST_TEST_MESSAGE("Testing angles (alpha, beta, gamma) = ("
			<< alpha << "," << beta << "," << gamma << ")");
	test_wigner_D_angles(alpha, beta, gamma);
}

BOOST_AUTO_TEST_SUITE_END()
