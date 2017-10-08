/*	This file test_ReadVASPPoscar.cpp is part of elephon.
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
 *  Created on: May 14, 2017
 *      Author: A. Linscheid
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Input_test
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "IOMethods/ReadVASPPoscar.h"
#include "fixtures/MockStartup.h"
#include <assert.h>
#include <vector>

BOOST_AUTO_TEST_CASE( Read_postcar )
{
	elephon::test::fixtures::MockStartup ms;
	elephon::IOMethods::ReadVASPPoscar filerreader;
	std::vector<std::pair<std::string, double>> atomOrder(1, std::make_pair(std::string("Al"), 26.981));
	filerreader.read_file( (ms.get_data_for_testing_dir() / "Al" / "vasp" / "conventional" / "POSCAR").string(), atomOrder );

	BOOST_REQUIRE( filerreader.get_atoms_list().size() == 4 );

	for ( auto a : filerreader.get_atoms_list() )
	{
		BOOST_REQUIRE( a.get_kind().compare("Al") == 0 );
	}

	auto A = filerreader.get_lattice_matrix();
	BOOST_REQUIRE( A.size() == 9 );
	for ( int i = 0 ; i < 3; ++i)
		for ( int j = 0 ; j < 3; ++j)
			BOOST_REQUIRE( std::fabs( A[i*3+j] - (i==j?4.038930:0.0) ) < 1e-6 );

}

BOOST_AUTO_TEST_CASE( Read_Al_fcc_primitive_vasp )
{
	elephon::test::fixtures::MockStartup ms;
	elephon::IOMethods::ReadVASPPoscar filerreader;
	std::vector<std::pair<std::string, double>> atomOrder(1, std::make_pair(std::string("Al"), 26.981));
	filerreader.read_file( (ms.get_data_for_testing_dir() / "Al" / "vasp" / "fcc_primitive" / "POSCAR").string(),
			atomOrder);

	BOOST_REQUIRE( filerreader.get_atoms_list().size() == 1 );

	for ( auto a : filerreader.get_atoms_list() )
	{
		BOOST_CHECK( a.get_kind().compare("Al") == 0 );
	}

	std::vector<double> latmatRef = { 	2.8560000000, 1.4280000000, 1.4280000000,
										0.0000000000, 2.4733685532, 0.8244561844,
										0.0000000000, 0.0000000000, 2.3319142351};

	auto A = filerreader.get_lattice_matrix();
	BOOST_REQUIRE( A.size() == 9 );
	for ( int i = 0 ; i < 3; ++i)
		for ( int j = 0 ; j < 3; ++j)
			BOOST_CHECK_CLOSE( A[i*3+j] , latmatRef[i*3+j], 0.0000001 );
}
