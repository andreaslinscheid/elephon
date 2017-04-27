/*	This file test_Input.cpp is part of elephon.
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
 *  Created on: Apr 25, 2017
 *      Author: A. Linscheid
 */
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Input_test
#include <boost/test/unit_test.hpp>
#include <string>
#include "IOMethods/Input.h"

BOOST_AUTO_TEST_CASE( Default_Args )
{
	int argc = 1;
	char const * argv[argc+1];
	argv[0] = "unused binary name";
	argv[1] = "../IOMethods/test_input_file.dat";
	elephon::IOMethods::Input input(argc,argv);
	auto sc = input.get_super_cell_size();
	BOOST_CHECK( sc.size() == 3 );
	BOOST_CHECK( sc[0] == 2 );
	BOOST_CHECK( sc[1] == 2 );
	BOOST_CHECK( sc[2] == 1 );

	BOOST_CHECK( input.get_numFS() == 2000 );
}
