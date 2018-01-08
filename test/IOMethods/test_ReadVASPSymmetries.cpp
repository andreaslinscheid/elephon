/*	This file test_ReadVASPSymmetries.cpp is part of elephon.
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
#include "IOMethods/ReadVASPSymmetries.h"
#include "fixtures/MockStartup.h"
#include <string>
#include <vector>

BOOST_AUTO_TEST_SUITE( ReadVASPSymmetries )

BOOST_AUTO_TEST_CASE( Read_symmetries )
{
	elephon::test::fixtures::MockStartup ms;
	auto testd = ms.get_data_for_testing_dir() / "Al" / "vasp" / "conventional";
	elephon::IOMethods::ReadVASPSymmetries symReader;
	symReader.read_file( (testd / "OUTCAR").string() );

	BOOST_REQUIRE( symReader.get_symmetries().size() == 48*9 );

	int refIndex  = 22;
	std::vector<int> ref_sym_22 = {0,0,1,0,-1,0,1,0,0};
	for ( int i = 0; i < 3 ; ++i)
		for ( int j = 0; j < 3 ; ++j)
			BOOST_REQUIRE(symReader.get_symmetries()[(refIndex*3+i)*3+j] == ref_sym_22[i*3+j]);

	BOOST_REQUIRE( symReader.get_fractionTranslations().size() == 48*3 );
}

BOOST_AUTO_TEST_SUITE_END()
