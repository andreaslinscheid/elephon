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
#include "fixtures/MockStartup.h"

BOOST_AUTO_TEST_CASE( Print_Manual )
{
	elephon::test::fixtures::MockStartup ms;
	auto mainDir = ms.get_data_for_testing_dir() / ".." / ".." ;
	auto manualPath = mainDir / "manual.txt";

	elephon::IOMethods::InputOptions options;
	options.build_input_manual( manualPath.string() );
}

BOOST_AUTO_TEST_CASE( Default_Args )
{
	using namespace boost::filesystem;

	elephon::test::fixtures::MockStartup ms;
	auto testd = ms.get_data_for_testing_dir() / "phony" ;

	elephon::IOMethods::InputOptions options;
	ms.simulate_elephon_input( (testd / "test_input_file.dat").string(), "scell=2 2 1\nnumFS=2000", options);

	BOOST_CHECK( options.get_scell().size() == 3 );
	BOOST_CHECK( options.get_scell()[0] == 2 );
	BOOST_CHECK( options.get_scell()[1] == 2 );
	BOOST_CHECK( options.get_scell()[2] == 1 );

	BOOST_CHECK( options.get_numFS() == 2000 );
}
