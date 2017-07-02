/*	This file test_LinearAlgebraInterface.cpp is part of elephon.
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
 *  Created on: Jun 16, 2017
 *      Author: A. Linscheid
 */
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Input_test
#include <boost/test/unit_test.hpp>
#include "Algorithms/LinearAlgebraInterface.h"
#include <iostream>
#include <vector>
#include <complex>

void print_mat( std::vector<double> const & v, int n , int m)
{
	for ( int i = 0 ; i < n ; ++i )
	{
		for ( int j = 0 ; j < m ; ++j )
			std::cout << v[i*n+j] << '\t';
		std::cout << '\n';
	}
	std::cout << std::endl;
}

BOOST_AUTO_TEST_CASE( Diagonalization )
{
	//Reproduce numpy
	auto ii = std::complex<double>(0,1);
	std::complex<double> sq2 = 1.0/std::sqrt(2);
	std::complex<double> sq2i = ii/std::sqrt(2);
	std::vector< std::complex<double> > A = { sq2, sq2i, -sq2i, -sq2 };

	std::vector<double> refEVal = { -1,   1 };
	std::vector< std::complex<double> > refEVec = { 1 , -ii , -ii , 1 };

	elephon::Algorithms::LinearAlgebraInterface linalg;
	std::vector<double> resultEVal;
	std::vector< std::complex<double> > resultEVec;
	linalg.diagonalize_hermitian( true, true, A ,resultEVec, resultEVal);

	std::sort( resultEVal.begin(), resultEVal.end() );
	BOOST_REQUIRE( resultEVec.size() == refEVec.size() );
	BOOST_REQUIRE( resultEVal.size() == refEVal.size() );
	for ( int i = 0 ; i < resultEVal.size(); ++i )
		BOOST_CHECK_CLOSE( refEVal[i], resultEVal[i], 0.000001 );
}

BOOST_AUTO_TEST_CASE( Pseudo_Inverse )
{
	//Reproduce numpy
	std::vector<double> A = { 1,2,3,4,5,6, 7,8,9,10,11,24, 13,14,15,16,17,18 };

	std::vector<double> ref = { -2.83333333e-01,   3.33333333e-02,   5.00000000e-02,
							    -1.45833333e-01,   8.33333333e-03,   3.75000000e-02,
							    -8.33333333e-03,  -1.66666667e-02,   2.50000000e-02,
							     1.29166667e-01,  -4.16666667e-02,   1.25000000e-02,
							     2.66666667e-01,  -6.66666667e-02,  -5.89805982e-17,
							    -4.16666667e-02,   8.33333333e-02,  -4.16666667e-02 };

	elephon::Algorithms::LinearAlgebraInterface linalg;
	std::vector<double> result;
	linalg.pseudo_inverse( A , 3, 6, result );

	BOOST_REQUIRE( result.size() == ref.size() );
	BOOST_REQUIRE( result.size() == A.size() );
	for ( int i = 0 ; i < A.size(); ++i )
		BOOST_CHECK_SMALL( result[i]-ref[i], 1e-6 );
}

// 48*2 x 3 matrix from a real world example
void real_world_data(std::vector<double> & A, int & m, int & n ){
	A = std::vector<double>({1,0,0,-1,0,-0,-1,0,0,1,0,0,0.5,0.288675,0.816497,-0.5,-0.288675,-0.816497,
			-0.5,-0.288675,-0.816497,0.5,0.288675,0.816497,0.5,0.866025,0,-0.5,-0.866025,0,-0.5,-0.866025,0,0.5,0.866025,
			0,0.816497,0.471405,0.333333,-0.816497,-0.471405,-0.333333,-0.816497,-0.471405,-0.333333,0.816497,0.471405,
			0.333333,0.816497,0.471405,0.333333,-0.816497,-0.471405,-0.333333,-0.816497,-0.471405,-0.333333,0.816497,
			0.471405,0.333333,0.816497,0.471405,0.333333,-0.816497,-0.471405,-0.333333,-0.816497,-0.471405,-0.333333,
			0.816497,0.471405,0.333333,0.5,0.866025,0,-0.5,-0.866025,0,-0.5,-0.866025,0,0.5,0.866025,-0,1,0,0,-1,0,0,
			-1,0,0,1,0,0,0.5,0.288675,0.816497,-0.5,-0.288675,-0.816497,-0.5,-0.288675,-0.816497,0.5,0.288675,0.816497,
			-0.5,-0.288675,-0.816497,0.5,0.288675,0.816497,0.5,0.288675,0.816497,-0.5,-0.288675,-0.816497,-0.5,-0.866025,
			0,0.5,0.866025,0,0.5,0.866025,0,-0.5,-0.866025,0,-1,0,0,1,0,0,1,0,0,-1,0,0,-0.816497,-0.471405,-0.333333,
			0.816497,0.471405,0.333333,0.816497,0.471405,0.333333,-0.816497,-0.471405,-0.333333,-0.816497,-0.471405,
			-0.333333,0.816497,0.471405,0.333333,0.816497,0.471405,0.333333,-0.816497,-0.471405,-0.333333,-0.816497,
			-0.471405,-0.333333,0.816497,0.471405,0.333333,0.816497,0.471405,0.333333,-0.816497,-0.471405,-0.333333,
			-0.5,-0.866025,0,0.5,0.866025,0,0.5,0.866025,0,-0.5,-0.866025,0,-1,0,0,1,0,-0,1,0,0,-1,0,0,-0.5,-0.288675,
			-0.816497,0.5,0.288675,0.816497,0.5,0.288675,0.816497,-0.5,-0.288675,-0.816497,1,0,0,-1,0,0,-1,0,0,1,0,0,
			0.5,0.288675,0.816497,-0.5,-0.288675,-0.816497,-0.5,-0.288675,-0.816497,0.5,0.288675,0.816497,0.5,0.866025,
			0,-0.5,-0.866025,0,-0.5,-0.866025,0,0.5,0.866025,0,-1,0,0,1,0,0,1,0,0,-1,-0,0,-0.5,-0.866025,0,0.5,0.866025,
			0,0.5,0.866025,0,-0.5,-0.866025,-0,-0.5,-0.288675,-0.816497,0.5,0.288675,0.816497,0.5,0.288675,0.816497,
			-0.5,-0.288675,-0.816497});
	m = 48*2;
	n = 3;
}

BOOST_AUTO_TEST_CASE( SVD )
{
	int m,n;
	std::vector<double>A;
	real_world_data(A,m,n);

	std::vector<double> ref_sv = {8.4852824 ,  3.46410317,  3.46410038};

	elephon::Algorithms::LinearAlgebraInterface linalg;
	std::vector<double> result, U,VT;
	assert( A.size() == m*n );
	linalg.svd( A , m, n, U,VT,result);

	BOOST_CHECK_EQUAL( result.size() , 3 );
	BOOST_REQUIRE( result.size() == ref_sv.size() );
	for ( int i = 0 ; i < result.size(); ++i )
		BOOST_CHECK_SMALL( result[i]-ref_sv[i], 1e-6 );
}

BOOST_AUTO_TEST_CASE( Nullspace )
{
	std::vector<double> A = { 1,2,3,4,5,6, 7,8,9,10,11,12, 13,14,15,16,17,18 };
	std::vector<double> AT = A;
	for ( int i = 0 ; i < 3; ++i)
		for ( int j = 0 ; j < 6; ++j)
			AT[j*3+i] = A[i*6+j];

	std::vector<double> ref = { 0.40824829, -0.81649658, 0.40824829};

	elephon::Algorithms::LinearAlgebraInterface linalg;
	std::vector<double> result;
	int kdim;
	linalg.null_space( A , 6, 3, kdim, result , 1e-5);

	BOOST_CHECK_EQUAL( kdim , 1 );
	BOOST_REQUIRE( result.size() == ref.size() );
	for ( int i = 0 ; i < result.size(); ++i )
		BOOST_CHECK_SMALL( result[i]-ref[i], 1e-6 );

	int m,n;
	real_world_data(A,m,n);
	linalg.null_space( A , m, n, kdim, result , 1e-5);
	BOOST_CHECK_EQUAL( kdim , 0 );
}