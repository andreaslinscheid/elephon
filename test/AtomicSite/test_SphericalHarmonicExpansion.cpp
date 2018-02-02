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
#include "AtomicSite/EulerAngles.h"
#include "AtomicSite/WignerDMatrix.h"
#include "AtomicSite/SphericalHarmonicExpansion.h"
#include "Auxillary/AlignedVector.h"
#include "Auxillary/memory_layout_functions.hpp"
#include "AtomicSite/RadialGrid.h"
#include "LatticeStructure/RegularBareGrid.h"
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <vector>
#include <complex>

BOOST_AUTO_TEST_SUITE( SphericalHarmonicExpansion )

elephon::AtomicSite::SphericalHarmonicExpansion
create_radial_constant_data(
		std::vector<int> ls,
		std::vector<int> ms,
		int lMax,
		int RMax)
{
	elephon::AtomicSite::RadialGrid rgrid;
	std::vector<double> points(RMax);
	double radius = 2.0*M_PI/std::pow(3.0,1.0/3.0);
	for (int ip = 0 ; ip < RMax ; ++ip)
		points[ip] = radius*static_cast<double>(ip+1)/static_cast<double>(RMax);
	rgrid.initialize({0.0, 0.0, 0.0}, radius, std::move(points));
	int nElem = (lMax+1)*(lMax+1)*RMax;
	elephon::Auxillary::alignedvector::ZV constant_data(nElem, 0.0);
	elephon::AtomicSite::SphericalHarmonicExpansion test_data;

	// now, only fill the radial data of constant '1' for all l and m.
	for (auto l : ls)
		for (auto m: ms)
			std::fill(constant_data.begin() + elephon::Auxillary::memlayout::angular_momentum_layout(l,m)*RMax,
					constant_data.begin() +  elephon::Auxillary::memlayout::angular_momentum_layout(l,m)*RMax + RMax,
					1.0);

	test_data.initialize(lMax, std::move(constant_data), std::move(rgrid));
	return test_data;
}

// This function checks if a SphericalHarmonicExpansion matches
// a cosine(2pi z) function in 3D of a 5^3 grid up to numerical precision of 1e-5
template<class functor>
void
check_is_function_on_grid_125(
		elephon::AtomicSite::SphericalHarmonicExpansion const & she,
		functor const & f)
{
	elephon::LatticeStructure::RegularBareGrid testGrid;
	testGrid.initialize({5, 5, 5});
	auto gridVectors = testGrid.get_all_vectors_grid();
	elephon::Auxillary::alignedvector::ZV interpolatedData;
	she.interpolate(gridVectors, interpolatedData);
	assert(interpolatedData.size() == testGrid.get_num_points());
	double diff = 0.0;
	for (int iGP = 0 ; iGP < testGrid.get_num_points(); ++iGP)
	{
		diff += std::abs(f(testGrid.get_vector_direct(iGP))-interpolatedData[iGP])
				/static_cast<double>(testGrid.get_num_points());
	}
	BOOST_CHECK_SMALL(diff, 1e-5);
}

BOOST_AUTO_TEST_CASE( test_she_datalayout )
{
	// first check if the data layout is correct.
	elephon::AtomicSite::RadialGrid rgrid;
	elephon::AtomicSite::SphericalHarmonicExpansion layout_test;
	elephon::Auxillary::alignedvector::ZV compareData;
	for (int l = 0; l<=5 ; ++l)
		for (int m = -l; m <= l; ++m)
			compareData.push_back(std::complex<double>(l,m));
	layout_test.initialize(5, compareData, rgrid);

	for (int l = 0; l<=5 ; ++l)
		for (int m = -l; m <= l; ++m)
			BOOST_CHECK_SMALL(std::abs(layout_test(0,m,l)-std::complex<double>(l,m)), 1e-10);
}

BOOST_AUTO_TEST_CASE( test_cosine_function )
{
	// now check that we recover a simple cos(2pi*x) function if we specify an expansion
	// that has only a component in Y_l=1,m=0 channel and a constant radial component.
	auto cos_z_test = create_radial_constant_data({1}, {0}, /*l_max = */5, /*Rmax = */50);

	// The above is fancy way of representing a cos(2pi*z) function in 3D. Check that indeed this is what we get.
	auto ylm_10 = [] (std::vector<double> const & v)
		{
			double z = v[2];
			double r = std::sqrt(std::pow(v[0],2) + std::pow(v[1],2) + std::pow(v[2],2));
			if ( r < 1e-10 )
				return std::complex<double>(std::sqrt(3.0/4.0/M_PI));
			return std::complex<double>(z/r)*std::sqrt(3.0/4.0/M_PI);
		};
	check_is_function_on_grid_125(cos_z_test, ylm_10);

	// now also check l=1 m=-1
	auto ylm_1m1 = create_radial_constant_data({1}, {-1}, /*l_max = */5, /*Rmax = */50);
	check_is_function_on_grid_125(ylm_1m1,
			[](std::vector<double> const & v)
			{
				double r = std::sqrt(std::pow(v[0],2) + std::pow(v[1],2) + std::pow(v[2],2));
				if (r < 1e-10)
					return std::complex<double>(0);
				return std::complex<double>(v[0]/r, -v[1]/r)*std::sqrt(3.0/8.0/M_PI);
			});

	// now also check l=1 m=+1
	auto ylm_1p1 = create_radial_constant_data({1}, {+1}, /*l_max = */5, /*Rmax = */50);
	check_is_function_on_grid_125(ylm_1p1,
			[](std::vector<double> const & v)
			{
				double r = std::sqrt(std::pow(v[0],2) + std::pow(v[1],2) + std::pow(v[2],2));
				if (r < 1e-10)
					return std::complex<double>(0);
				return std::complex<double>(-v[0]/r, -v[1]/r)*std::sqrt(3.0/8.0/M_PI);
			});

	// Now check a linear combination
	auto cos_x_y_test = create_radial_constant_data({1}, {1,-1}, /*l_max = */5, /*Rmax = */50);

	// Check that the resulting function is Y_1,-1+Y_1,1
	check_is_function_on_grid_125(cos_x_y_test, [] (std::vector<double> const & v)
			{
				double r = std::sqrt(std::pow(v[0],2) + std::pow(v[1],2) + std::pow(v[2],2));
				if ( r < 1e-10 )
					return std::complex<double>(0);
				return std::complex<double>(0, -v[1]/r)*std::sqrt(3.0/2.0/M_PI);
			} );
}

BOOST_AUTO_TEST_CASE( test_rotation )
{
	// Produce an initial ylm_1,-1 and rotate it.
	// rotating around the y axis by 180Deg has to transform the result into -ylm_1,1
	// This tests: 	1. Decomposition of a rotation matrix into Euler angles.
	//				2. Construction of the Wigner D rotation matrix.
	//				3. Application of a Wigner D rotation matrix on the spherical harmoic
	//					expansion. Step 2 is not verified since it is are subject
	//					to independent unit tests.
	std::vector<double> rotationMatrix{-1.0, 0.0, 0.0,
										0.0, 1.0, 0.0,
									    0.0, 0.0,-1.0 };

	double alpha, beta, gamma;
	elephon::AtomicSite::eulerAngles(rotationMatrix, alpha, beta, gamma);

	BOOST_CHECK_SMALL(alpha, 1e-5);
	BOOST_CHECK_SMALL(beta-M_PI, 1e-5);
	BOOST_CHECK_SMALL(gamma, 1e-5);

	std::vector<elephon::AtomicSite::WignerDMatrix> wd(6);
	for (int l = 0; l <= 5; ++l)
		wd[l].initialize(l, alpha, beta, gamma);
	auto rotOp = std::make_shared<decltype(wd)>(std::move(wd));

	// construct a formal symmetry operation. In this test, only the angular part is used.
	std::vector<int> rotMatLatticeBasis{-1, 0, 0, 0, 1, 0, 0, 0, -1};
	std::vector<double> fracTransZero{0.0, 0.0, 0.0};
	elephon::symmetry::SymmetryOperation sop(rotMatLatticeBasis.begin(), rotMatLatticeBasis.end(),
											 fracTransZero.begin(), fracTransZero.end(),
											 rotationMatrix.begin(), rotationMatrix.end(),
											 fracTransZero.begin(), fracTransZero.end(),
											 rotOp);

	// Create l=1 m=-1
	auto ylm_1m1 = create_radial_constant_data({1}, {-1}, /*l_max = */5, /*Rmax = */50);

	// Rotate by 180 Deg about the y axis
	ylm_1m1.transform(sop);

	// now check that this equals -Ylm for l=1 m=+1
	check_is_function_on_grid_125(ylm_1m1,
			[](std::vector<double> const & v)
			{
				double r = std::sqrt(std::pow(v[0],2) + std::pow(v[1],2) + std::pow(v[2],2));
				if (r < 1e-10)
					return std::complex<double>(0);
				return std::complex<double>(-v[0]/r, -v[1]/r)*std::sqrt(3.0/8.0/M_PI);
			});
}

BOOST_AUTO_TEST_CASE( test_data_fit )
{
	// in  this test, we generate data that is a known function of radius
	// times a spherical harmonic. We then confirm that with a fit, the coefficients are matching
	// this behavior.
	struct testDataGenerator
	{
		testDataGenerator(int l, int m, double a, std::vector<double> center) : l_(l), m_(m), a_(a), c_(std::move(center)) {
			assert(c_.size()==3);
		};

		void interpolate( std::vector<double> & coordinates, std::complex<double> * data) const
		{
			const int np = coordinates.size()/3;
			for (int ip = 0 ; ip < np; ++ip )
			{
				double x = coordinates[ip*3+0] - c_[0];
				double y = coordinates[ip*3+1] - c_[1];
				double z = coordinates[ip*3+2] - c_[2];
				double r, theta, phi;
				elephon::Algorithms::helperfunctions::compute_spherical_coords(x, y, z, r, theta, phi);
				data[ip] = boost::math::spherical_harmonic(l_, m_, theta, phi )*(1.0 - a_*r);
			}
		}

		int l_, m_;
		double a_;
		std::vector<double> c_;
	};

	const int RMax = 50;
	const int LMax = 10;
	const int LTest = 1;
	const int MTest = 1;
	const double aTest = 0.5;

	std::vector<double> center{0.5, 0.25, 0.125};
	testDataGenerator gen(LTest, MTest, aTest, center);

	elephon::AtomicSite::RadialGrid rgrid;
	std::vector<double> points(RMax);
	double radius = 2.0*M_PI/std::pow(3.0,1.0/3.0);
	for (int ip = 0 ; ip < RMax ; ++ip)
		points[ip] = radius*static_cast<double>(ip+1)/static_cast<double>(RMax);
	rgrid.initialize(center, radius, std::move(points));

	elephon::AtomicSite::SphericalHarmonicExpansion shexp;
	shexp.fit_to_data(gen, LMax, rgrid);

	for (int iL = 0 ; iL < LMax; ++iL)
		for (int iM = -iL ; iM <= iL; ++iM)
		{
			double diff = 0.0;
			for (int iR = 0 ; iR < RMax; ++iR)
			{
				double r = rgrid.get_radius(iR);
				std::complex<double> expected = ((iL == LTest)&&(iM == MTest) ?
											std::complex<double>(1.0 - aTest*r) : std::complex<double>(0.0));
				std::complex<double> obtained = shexp(iR, iM, iL);
				diff += std::abs(obtained - expected)/RMax;
			}
			BOOST_CHECK_SMALL(diff, 1e-6);
		}
}

BOOST_AUTO_TEST_SUITE_END()
