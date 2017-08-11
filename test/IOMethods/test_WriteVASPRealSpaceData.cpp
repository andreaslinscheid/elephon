/*	This file test_WriteVASPRealSpaceData.cpp is part of elephon.
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
 *  Created on: Jul 3, 2017
 *      Author: A. Linscheid
 */


#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE output_VASP
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "IOMethods/WriteVASPRealSpaceData.h"
#include "IOMethods/ReadVASPLocpot.h"
#include "fixtures/MockStartup.h"
#include "fixtures/DataLoader.h"
#include <string>
#include <vector>

BOOST_AUTO_TEST_CASE( MgB2_realspace_data_format )
{
	test::fixtures::MockStartup ms;
	auto testf = ms.get_data_for_testing_dir() / "vasp_realspace_data.dat";

	auto loadD = ms.get_data_for_testing_dir() / "MgB2" / "vasp" / "ldos";
	std::string content = std::string("root_dir=")+loadD.string()+"\n";
	test::fixtures::DataLoader dl;
	auto unitCell = dl.load_unit_cell( content );

	elephon::IOMethods::WriteVASPRealSpaceData wd;
	std::vector<int> dims{20,20,24};
	std::vector<double> data(dims[0]*dims[1]*dims[2], 0.0);

	wd.write_file( testf.string(), "Test file comment", dims, unitCell, data );

	elephon::IOMethods::ReadVASPLocpot locpot;
	locpot.read_scf_potential(testf.string(), dims, data);
}

BOOST_AUTO_TEST_CASE( write_VASP_realspace_output )
{
	test::fixtures::MockStartup ms;
	auto testf = ms.get_data_for_testing_dir() / "vasp_realspace_data.dat";

	auto loadD = ms.get_data_for_testing_dir() / "Al" / "vasp" / "fcc_primitive";
	std::string content = std::string("root_dir=")+loadD.string()+"\n";
	test::fixtures::DataLoader dl;
	auto unitCell = dl.load_unit_cell( content );

	elephon::IOMethods::WriteVASPRealSpaceData wd;
	std::vector<int> dims{6,6,6};
	std::vector<double> data(dims[0]*dims[1]*dims[2]);
	for ( int i = 0 ; i < dims[0]*dims[1]*dims[2] ; ++i )
		data[i] = i-2;

	//Also put some corner cases for the formating into the file
	data[5] = 1e99;
	data[5] = -1e-100;
	wd.write_file( testf.string(), "Test file comment", dims, unitCell, data );
}

