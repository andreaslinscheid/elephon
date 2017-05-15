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
#include <assert.h>
#include <vector>

BOOST_AUTO_TEST_CASE( Read_postcar )
{
	boost::filesystem::path p(__FILE__);
	boost::filesystem::path dir = p.parent_path();
	std::string Al_test_poscar = std::string(dir.c_str())+"/../IOMethods/POSCAR_Al_test.dat";

	elephon::IOMethods::ReadVASPPoscar filerreader;
	filerreader.read_file( Al_test_poscar );

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
