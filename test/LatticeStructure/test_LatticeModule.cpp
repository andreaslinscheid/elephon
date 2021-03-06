/*	This file test_LatticeModule.cpp is part of elephon.
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
 *  Created on: May 16, 2017
 *      Author: A. Linscheid
 */
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "LatticeStructure/LatticeModule.h"

BOOST_AUTO_TEST_SUITE( LatticeModule )

BOOST_AUTO_TEST_CASE( LatticeModule_BasisRelation )
{
	std::vector<double> latticeMatrix = {	4.038930, 0.000000, 0.000000,
											0.000000, 4.038930, 0.000000,
											0.000000, 0.000000, 4.038930};

	elephon::LatticeStructure::LatticeModule lattice;
	lattice.initialize(latticeMatrix);

	//Check that the
	std::vector<double> prod(9,0.0);
	auto A = lattice.get_latticeMatrix();
	auto B = lattice.get_reciprocal_latticeMatrix();
	for ( int i = 0 ; i < 3; ++i)
		for ( int j = 0 ; j < 3; ++j)
			for ( int k = 0 ; k < 3; ++k)
				prod[i*3+j] += B[k*3+i]*A[k*3+j];

	for ( int i = 0 ; i < 3; ++i)
		for ( int j = 0 ; j < 3; ++j)
			BOOST_REQUIRE( std::fabs(prod[i*3+j] - (i==j?1.0:0)) < 1e-6 );
}
BOOST_AUTO_TEST_SUITE_END()
