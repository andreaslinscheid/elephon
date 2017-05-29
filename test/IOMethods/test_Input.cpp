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
#include <boost/filesystem.hpp>
#include <string>
#include "IOMethods/Input.h"

BOOST_AUTO_TEST_CASE( Default_Args )
{
	using namespace boost::filesystem;

	char * prog = strdup("program name");
	char * arg = strdup( (path(__FILE__).parent_path() / "test_input_file.dat" ).c_str());
	char *argv[] = {prog, arg, NULL};
	int argc = sizeof(argv) / sizeof(char*) - 1;
	elephon::IOMethods::Input input(argc,argv);
	elephon::IOMethods::InputOptions options = input.get_opts();
	delete [] prog;
	delete [] arg;


	BOOST_CHECK( options.get_scell().size() == 3 );
	BOOST_CHECK( options.get_scell()[0] == 2 );
	BOOST_CHECK( options.get_scell()[1] == 2 );
	BOOST_CHECK( options.get_scell()[2] == 1 );

	BOOST_CHECK( options.get_numFS() == 2000 );
}
