/*	This file test_AtomicSite.cpp is part of elephon.
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
 *  Created on: Jan 2, 2018
 *      Author: A. Linscheid
 */

#define BOOST_TEST_MODULE Input_test
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "fixtures/MockStartup.h"
#include "AtomicSite/EulerAngles.h"
#include <vector>

BOOST_AUTO_TEST_SUITE( ASSymmetry_test )

BOOST_AUTO_TEST_CASE( test_euler_angle_decomp )
{
	// Constructed Euler matrix for (alpha, beta, gamma) = (Pi/2, Pi/4, Pi)
	double s2 = 1.0/std::sqrt(2.0);
	std::vector<double> EulerMatrix = {  0.0, 1.0, 0.0,
										-s2,  0.0, s2,
										 s2,  0.0, s2};

	double alpha, beta, gamma;
	elephon::AtomicSite::eulerAngles(EulerMatrix, alpha, beta, gamma);

	BOOST_CHECK_CLOSE(alpha, M_PI/2.0, 1e-10);
	BOOST_CHECK_CLOSE(beta, M_PI/4.0, 1e-10);
	BOOST_CHECK_CLOSE(gamma, M_PI, 1e-10);

	// Corner cases beta = 0: (alpha, beta, gamma) = (Pi/2, 0, Pi/3)
	//  NOTE: Since beta is 0, this is effectively a rotation twice about the z axis.
	//  This has not a unique solution, so by convention only alpha is rotated as pi/3+pi/2=5/6 pi.
	EulerMatrix = std::vector<double>{ -std::sqrt(3.0)/2.0, -0.5, 				 0.0,
										0.5, 				-std::sqrt(3.0)/2.0, 0.0,
										0.0,  				0.0, 				 1.0};
	elephon::AtomicSite::eulerAngles(EulerMatrix, alpha, beta, gamma);
	BOOST_CHECK_CLOSE(alpha, 5.0/6.0*M_PI, 1e-10);
	BOOST_CHECK_CLOSE(beta, 0.0, 1e-10);
	BOOST_CHECK_CLOSE(gamma, 0.0, 1e-10);

	// Corner cases beta = Pi: (alpha, beta, gamma) = (Pi/2, Pi, Pi/3)
	//  NOTE: Since beta is Pi, this is effectively a rotation and then a rotation back around the z axis.
	//  This has not a unique solution, so by convention only alpha is rotated as pi/3-pi/2=-1/6 pi.
	EulerMatrix = std::vector<double>{ -std::sqrt(3.0)/2.0, 0.5, 				0.0,
										0.5, 				std::sqrt(3.0)/2.0, 0.0,
										0.0,  				0.0, 				-1.0};
	elephon::AtomicSite::eulerAngles(EulerMatrix, alpha, beta, gamma);
	BOOST_CHECK_CLOSE(alpha, -1.0/6.0*M_PI, 1e-10);
	BOOST_CHECK_CLOSE(beta, M_PI, 1e-10);
	BOOST_CHECK_CLOSE(gamma, 0.0, 1e-10);
}

BOOST_AUTO_TEST_SUITE_END()
