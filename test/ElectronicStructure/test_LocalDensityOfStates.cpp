/*	This file test_LocalDensityOfStates.cpp is part of elephon.
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
 *  Created on: Jul 26, 2017
 *      Author: A. Linscheid
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Input_test
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "IOMethods/VASPInterface.h"
#include "ElectronicStructure/LocalDensityOfStates.h"
#include "fixtures/MockStartup.h"
#include <vector>
#include <complex>
#include <cmath>

BOOST_AUTO_TEST_CASE( ldos_MgB2_vasp )
{
	test::fixtures::MockStartup ms;
	auto testd = ms.get_data_for_testing_dir() / "MgB2" / "vasp" / "ldos";
	auto outfile = testd / "ldos.dat";
	boost::filesystem::remove(outfile);

	std::string input = std::string()+
			"root_dir = "+testd.string()+"\n"
			"f_ldos = "+outfile.string()+"\n"
			"fftd = 64 64 64\n"
			"eldos = 0.0\n";
	elephon::IOMethods::InputOptions opts;
	ms.simulate_elephon_input(
			(testd / "infile").string(),
			input,
			opts);

	auto loader = std::make_shared<elephon::IOMethods::VASPInterface>(opts);

	elephon::ElectronicStructure::LocalDensityOfStates ldos;
	ldos.compute_ldos(opts.get_eldos(), loader);
	ldos.write_file(opts.get_f_ldos());

	BOOST_REQUIRE( boost::filesystem::exists(outfile) );
	boost::filesystem::remove(outfile);
}

BOOST_AUTO_TEST_CASE( ldos_Al_vasp )
{
	test::fixtures::MockStartup ms;
	auto testd = ms.get_data_for_testing_dir() / "Al" / "vasp" / "fcc_primitive";
	auto outfile = testd / "ldos.dat";
	boost::filesystem::remove(outfile);

	std::string input = std::string()+
			"root_dir = "+testd.string()+"\n"
			"f_ldos = "+outfile.string()+"\n"
			"fftd = 128 128 128\n"
			"eldos = 0.0\n";
	elephon::IOMethods::InputOptions opts;
	ms.simulate_elephon_input(
			(testd / "infile").string(),
			input,
			opts);

	auto loader = std::make_shared<elephon::IOMethods::VASPInterface>(opts);

	elephon::ElectronicStructure::LocalDensityOfStates ldos;
	ldos.compute_ldos(opts.get_eldos(), loader);
	ldos.write_file(opts.get_f_ldos());

	BOOST_REQUIRE( boost::filesystem::exists(outfile) );
	boost::filesystem::remove(outfile);
}
