/*	This file test_rotation_matrix_connecting_unit_vectors.cpp is part of elephon.
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
#include <boost/test/unit_test.hpp>
#include "Algorithms/rotation_matrix_connecting_unit_vectors.h"
#include <boost/multi_array.hpp>
#include <vector>
#include <cstdlib>
#include <ctime>

BOOST_AUTO_TEST_SUITE( rotation_matrix_connecting_unit_vectors )

std::vector<double>
get_random_finite_vector()
{
	std::vector<double> vec(3, 0.0);
	auto norm2 = [] (std::vector<double> const & a) {return a[0]*a[0]+a[1]*a[1]+a[2]*a[2];};
	while (norm2(vec) < 1e-8)
	{
		for (auto &xi : vec)
			xi = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
	}
	return vec;
}

std::vector<double>
get_random_unit_vector()
{
	auto vec = get_random_finite_vector();
	double norm = std::sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
	for (auto &v : vec)
		v /= norm;
	return vec;
}

BOOST_AUTO_TEST_CASE( identity )
{
	// verify that for identical vector the method yields the identity matrix
	srand (time(NULL));
	auto vec = get_random_unit_vector();

	BOOST_TEST_MESSAGE("Got random vector ("+std::to_string(vec[0])
				+","+std::to_string(vec[1])+","+std::to_string(vec[2])+")");

	boost::multi_array<double,2> rotationMatrix;
	elephon::Algorithms::rotation_matrix_connecting_unit_vectors(
			vec,
			vec,
			rotationMatrix);

	const double accu = 1e-8;
	BOOST_CHECK_SMALL(rotationMatrix[0][0]-1.0, accu);
	BOOST_CHECK_SMALL(rotationMatrix[0][1]-0.0, accu);
	BOOST_CHECK_SMALL(rotationMatrix[0][2]-0.0, accu);
	BOOST_CHECK_SMALL(rotationMatrix[1][0]-0.0, accu);
	BOOST_CHECK_SMALL(rotationMatrix[1][1]-1.0, accu);
	BOOST_CHECK_SMALL(rotationMatrix[1][2]-0.0, accu);
	BOOST_CHECK_SMALL(rotationMatrix[2][0]-0.0, accu);
	BOOST_CHECK_SMALL(rotationMatrix[2][1]-0.0, accu);
	BOOST_CHECK_SMALL(rotationMatrix[2][2]-1.0, accu);
}

BOOST_AUTO_TEST_CASE( rotate_pi )
{
	// test 2 special cases where v is rotated into -v. One case has x <= y and the other y>x
	std::vector<double> vec{0.1, 0.2, 0.3};
	double norm = std::sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
	for (auto &v: vec)
		v /= norm;
	auto minus_vec = vec;
	for (auto &ai : minus_vec)
		ai *= -1;
	boost::multi_array<double,2> rotationMatrix;
	elephon::Algorithms::rotation_matrix_connecting_unit_vectors(
			vec,
			minus_vec,
			rotationMatrix);

	// We obtain the matrix that generates the inverse vector from Mathematica, confirm that
	// this is the right matrix
	std::vector<double> reference_Matrix{	-1.0,  0.0, 				 0.0,
											 0.0,  0.3846153846153847, 	-0.9230769230769231,
											 0.0, -0.9230769230769231, 	-0.3846153846153847};
	auto check_ref_matrix = [] (std::vector<double> const & vec,
								std::vector<double> const & reference_Matrix) {
		for (int i = 0 ; i < 3; ++i)
		{
			double val = 0;
			for (int j = 0 ; j < 3; ++j)
				val += reference_Matrix[i*3+j]*vec[j];
			BOOST_CHECK_SMALL(val+vec[i], 1e-7);
		}
	};
	check_ref_matrix(vec, reference_Matrix);

	// Now check if that is the matrix we got from the method.
	for (int i = 0 ; i < 3; ++i)
		for (int j = 0 ; j < 3; ++j)
			BOOST_CHECK_SMALL(rotationMatrix[i][j]-reference_Matrix[i*3+j], 1e-8);

	// repeat for the other possible branch
	vec = std::vector<double>{0.2, 0.1, 0.3};
	norm = std::sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
	for (auto &v: vec)
		v /= norm;
	reference_Matrix = std::vector<double> { 0.3846153846153847,  0.0, -0.9230769230769231,
											 0.0, 				 -1.0,  0.0,
											-0.9230769230769231,  0.0, -0.3846153846153847};
	check_ref_matrix(vec, reference_Matrix);

	// main testing of the other branch
	minus_vec = vec;
	for (auto &ai : minus_vec)
		ai *= -1;
	elephon::Algorithms::rotation_matrix_connecting_unit_vectors(
			vec,
			minus_vec,
			rotationMatrix);

	// Again check if that is the matrix we got from the method.
	for (int i = 0 ; i < 3; ++i)
		for (int j = 0 ; j < 3; ++j)
			BOOST_CHECK_SMALL(rotationMatrix[i][j]-reference_Matrix[i*3+j], 1e-8);

}

// generate a random unit vector and rotate randomly around it
std::vector<double>
generate_random_rotation_matrix()
{
	auto phi = 2.0*M_PI*static_cast<double>(rand())/static_cast<double>(RAND_MAX);
	auto cos_phi = std::cos(phi);
	auto sin_phi = std::sin(phi);
	auto rot_vec = get_random_unit_vector();
	std::vector<double> R{
		cos_phi +std::pow(rot_vec[0],2)*(1-cos_phi),
		rot_vec[0]*rot_vec[1]*(1-cos_phi) - rot_vec[2]*sin_phi,
		rot_vec[0]*rot_vec[2]*(1-cos_phi) + rot_vec[1]*sin_phi,
		rot_vec[1]*rot_vec[0]*(1-cos_phi) + rot_vec[2]*sin_phi,
		cos_phi + std::pow(rot_vec[1],2)*(1-cos_phi),
		rot_vec[1]*rot_vec[2]*(1-cos_phi) - rot_vec[0]*sin_phi,
		rot_vec[2]*rot_vec[0]*(1-cos_phi) - rot_vec[1]*sin_phi,
		rot_vec[2]*rot_vec[1]*(1-cos_phi) + rot_vec[0]*sin_phi,
		cos_phi + std::pow(rot_vec[2],2)*(1-cos_phi)};
	return R;
}

void rot_check(std::vector<double> const & vec, std::vector<double> const & R)
{
	// apply matrix to generate rotated vector.
	std::vector<double> rot_vec(3, 0.0);
	for (int i = 0 ; i < 3 ; ++i)
	{
		for (int j = 0 ; j < 3 ; ++j)
			rot_vec[i] += R[i*3+j]*vec[j];
	}

	// pass vector and check that we get our rotation matrix back.
	boost::multi_array<double,2> rotationMatrix;
	elephon::Algorithms::rotation_matrix_connecting_unit_vectors(
			vec,
			rot_vec,
			rotationMatrix);

	// Again check if that is the matrix we got from the method.
	for (int i = 0 ; i < 3; ++i)
	{
		double val = 0.0;
		for (int j = 0 ; j < 3; ++j)
			val += rotationMatrix[i][j]*vec[j];
		BOOST_CHECK_SMALL(rot_vec[i]-val, 1e-8);
	}
}

BOOST_AUTO_TEST_CASE( rotate_general_example )
{
	std::vector<double> vec{0.538489,0.250337,0.804587};
	double norm = std::sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
	for (auto &v: vec)
		v /= norm;
	std::vector<double> R{-0.586623,0.498852,-0.637981,-0.322528,0.578694,0.749059, 0.742865,0.645182,-0.178581};
	rot_check(vec, R);
}

BOOST_AUTO_TEST_CASE( rotate_general )
{
	srand (time(NULL));
	auto vec = get_random_unit_vector();
	BOOST_TEST_MESSAGE("Got random vector ("+std::to_string(vec[0])
				+","+std::to_string(vec[1])+","+std::to_string(vec[2])+")");

	auto R = generate_random_rotation_matrix();
	BOOST_TEST_MESSAGE("Got random matrix \n("
				+std::to_string(R[0])+","+std::to_string(R[1])+","+std::to_string(R[2])+")\n"
				+std::to_string(R[3])+","+std::to_string(R[4])+","+std::to_string(R[5])+")\n"
				+std::to_string(R[6])+","+std::to_string(R[7])+","+std::to_string(R[8])+") ");
	rot_check(vec, R);
}

BOOST_AUTO_TEST_SUITE_END()
