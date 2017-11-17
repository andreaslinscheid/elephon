/*	This file test_KPath.cpp is part of elephon.
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
 *  Created on: Nov 1, 2017
 *      Author: A. Linscheid
 */

#define BOOST_TEST_MODULE Kpath_test
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "IOMethods/VASPInterface.h"
#include "IOMethods/ResourceHandler.h"
#include "fixtures/MockStartup.h"
#include "IOMethods/KPath.h"

BOOST_AUTO_TEST_CASE( KPath_simple )
{
	elephon::test::fixtures::MockStartup ms;
	auto testd = ms.get_data_for_testing_dir() / "phony";
	auto fkpath = testd / "kpath.dat";
	boost::filesystem::remove(fkpath);

	std::string input = std::string()+
			"root_dir = "+testd.string()+"\n"
			"f_kpath = "+fkpath.string()+"\n";
	elephon::IOMethods::InputOptions opts;
	ms.simulate_elephon_input(
			(testd / "infile").string(),
			input,
			opts);

	ms.write_kpath_file_fcc(fkpath.string());
	auto loader = std::make_shared<elephon::IOMethods::VASPInterface>(opts);
	auto resource = std::make_shared<elephon::IOMethods::ResourceHandler>(loader);
	auto kpath = resource->get_k_path();

	auto kpts = kpath->get_k_points();
	BOOST_CHECK_EQUAL( kpts.size()/3 , 40+40+40+57+1 );

	boost::filesystem::remove(fkpath);
}

