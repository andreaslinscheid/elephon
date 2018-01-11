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
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "fixtures/MockStartup.h"
#include "AtomicSite/EulerAngles.h"
#include "AtomicSite/WignerDMatrix.h"
#include "AtomicSite/SphericalHarmonicExpansion.h"
#include "Auxillary/AlignedVector.h"
#include "AtomicSite/RadialGrid.h"
#include "LatticeStructure/RegularBareGrid.h"
#include <vector>
#include <complex>
#include <iostream>
#include <cstdlib>
#include <ctime>

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

BOOST_AUTO_TEST_CASE( test_wigner_d_matrix )
{
	// This lambda tests the examples for l=0, l=1 and l=2 against
	// analytically known forms of the D matrix
	auto test_wigner_D_angles = [] (double alpha, double beta, double gamma)
	{

		elephon::AtomicSite::WignerDMatrix D0;
		D0.initialize(0, alpha, beta, gamma);
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
		D1.initialize(l, alpha, beta, gamma);
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
		D2.initialize(l, alpha, beta, gamma);
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

BOOST_AUTO_TEST_CASE( test_spherical_harmonic_expansion )
{
	// first check if the data layout is correct.
	elephon::AtomicSite::RadialGrid rgrid;
	elephon::AtomicSite::SphericalHarmonicExpansion layout_test;
	elephon::Auxillary::alignedvector::ZV compareData;
	for (int l = 0; l<=5 ; ++l)
		for (int m = -l; m <= l; ++m)
			compareData.push_back(std::complex<double>(l,m));
	layout_test.initialize(5, 1, compareData, rgrid);

	for (int l = 0; l<=5 ; ++l)
		for (int m = -l; m <= l; ++m)
			BOOST_CHECK_SMALL(std::abs(layout_test(0,m,l)-std::complex<double>(l,m)), 1e-10);

	// now check that we recover a simple cos(2pi*x) function if we specify an expansion
	// that has only a component in Y_l=1,m=0 channel and a constant radial component.
	std::vector<double> points(50);
	double radius = 2.0*M_PI/std::pow(3.0,1.0/3.0);
	for (int ip = 0 ; ip < 50 ; ++ip)
		points[ip] = radius*static_cast<double>(ip)/static_cast<double>(50);
	rgrid.initialize({0.0, 0.0, 0.0}, radius, std::move(points));
	int lmax = 5;
	int nElem = (lmax+1)*(lmax+1)*rgrid.get_num_R();
	elephon::Auxillary::alignedvector::ZV constant_data(nElem, 0.0);
	elephon::AtomicSite::SphericalHarmonicExpansion cos_x_test;
	// now, only fill the radial data of constant '1' for l=1 and m=0 with 1.
	std::fill(constant_data.begin() + cos_x_test.angular_momentum_layout(/*l =*/1, /*m =*/0),
			constant_data.begin() + cos_x_test.angular_momentum_layout(/*l =*/1, /*m =*/0) + points.size(),
			1.0);

	cos_x_test.initialize(5, 50.0, std::move(constant_data), std::move(rgrid));

	// The above is fancy way of representing a cos(2pi*z) function in 3D. Check that indeed this is what we get.
	elephon::LatticeStructure::RegularBareGrid testGrid;
	testGrid.initialize({50, 50, 50});
	auto gridVectors = testGrid.get_all_vectors_grid();
	elephon::Auxillary::alignedvector::ZV interpolatedData;
	cos_x_test.interpolate(gridVectors, interpolatedData);
	assert(interpolatedData.size() == testGrid.get_num_points());
	double diff = 0.0;
	for (int iGP = 0 ; iGP < testGrid.get_num_points(); ++iGP)
	{
		auto thisVector = testGrid.get_vector_direct(iGP);
		double z = thisVector[3];
		auto cos_z = std::sqrt(3.0/4.0/M_PI)*std::complex<double>(2.0*M_PI*z);
		diff += std::abs(cos_z-interpolatedData[iGP])
				/static_cast<double>(testGrid.get_num_points());
	}
	BOOST_CHECK_SMALL(diff, 1e-5);
}

BOOST_AUTO_TEST_SUITE_END()
