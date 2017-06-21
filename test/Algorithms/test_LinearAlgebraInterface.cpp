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

BOOST_AUTO_TEST_CASE( Pseudo_Inverse )
{
	//Reporduce numpy
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
		BOOST_CHECK_CLOSE( result[i], ref[i], 0.0001 );
}
