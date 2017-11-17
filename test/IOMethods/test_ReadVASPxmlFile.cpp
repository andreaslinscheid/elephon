/*	This file test_ReadVASPxmlFile.cpp is part of elephon.
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
 *  Created on: May 31, 2017
 *      Author: A. Linscheid
 */

#define BOOST_TEST_MODULE Input_test
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "IOMethods/ReadVASPxmlFile.h"
#include "fixtures/MockStartup.h"
#include <assert.h>
#include <vector>

BOOST_AUTO_TEST_CASE( Read_VASP_xml )
{
	elephon::test::fixtures::MockStartup ms;
	auto testd = ms.get_data_for_testing_dir() / "Al" / "vasp" / "conventional" / "supercell" ;
	elephon::IOMethods::ReadVASPxmlFile filerreader;
	filerreader.parse_file( (testd / "vasprun.xml").string() );

	double EFermi = 8.01896145;
	BOOST_REQUIRE( filerreader.get_Fermi_energy() == EFermi );

	std::vector<double> forces = {
		  0.03989949   ,    0.00000000   ,    0.00000000 ,
	     -0.00152643   ,    0.00353302   ,    0.00000000 ,
	      0.00161674   ,    0.00000000   ,    0.00000000 ,
	     -0.00266775   ,   -0.01081380   ,    0.00000000 ,
	     -0.01286095   ,    0.00000000   ,    0.00000000 ,
	     -0.00152643   ,   -0.00353302   ,    0.00000000 ,
	      0.00429365   ,    0.00000000   ,    0.00000000 ,
	     -0.00266775   ,    0.01081380   ,    0.00000000 ,
	     -0.01717443   ,    0.00000000   ,    0.00000000 ,
	     -0.00270537   ,   -0.00386870   ,    0.00000000 ,
	     -0.00589763   ,    0.00000000   ,    0.00000000 ,
	     -0.00253047   ,    0.01055867   ,    0.00000000 ,
	     -0.00310275   ,    0.00000000   ,    0.00000000 ,
	     -0.00270537   ,    0.00386870   ,    0.00000000 ,
	      0.01208592   ,    0.00000000   ,    0.00000000 ,
	     -0.00253047   ,   -0.01055867   ,    0.00000000 };

	BOOST_REQUIRE(filerreader.get_forces() == forces);
}

BOOST_AUTO_TEST_CASE( Read_VASP_xml_kpoints )
{
	elephon::test::fixtures::MockStartup ms;
	auto testd = ms.get_data_for_testing_dir() / "Al" / "vasp" / "conventional" / "supercell" ;
	elephon::IOMethods::ReadVASPxmlFile filerreader;
	filerreader.parse_file( (testd / "vasprun.xml").string() );

	std::vector<double> kpoints = {
			  0.00000000 , 0.00000000 , 0.08333333 ,
			  0.33333333 , 0.00000000 , 0.08333333 ,
			  0.00000000 , 0.33333333 , 0.08333333 ,
			  0.33333333 , 0.33333333 , 0.08333333 ,
			  0.00000000 , 0.00000000 , 0.25000000 ,
			  0.33333333 , 0.00000000 , 0.25000000 ,
			  0.00000000 , 0.33333333 , 0.25000000 ,
			  0.33333333 , 0.33333333 , 0.25000000 ,
			  0.00000000 , 0.00000000 , 0.41666667 ,
			  0.33333333 , 0.00000000 , 0.41666667 ,
			  0.00000000 , 0.33333333 , 0.41666667 ,
			  0.33333333 , 0.33333333 , 0.41666667 };

	BOOST_REQUIRE(filerreader.get_k_points() == kpoints);

}

BOOST_AUTO_TEST_CASE( Read_VASP_xml_kpoints_fcc_Al )
{
	elephon::test::fixtures::MockStartup ms;
	auto rootDir = ms.get_data_for_testing_dir() / "Al" / "vasp" / "fcc_primitive";
	elephon::IOMethods::ReadVASPxmlFile filerreader;
	filerreader.parse_file( (rootDir / "vasprun.xml").string() );

	auto p = std::minmax_element( filerreader.get_k_points().begin(), filerreader.get_k_points().end() );
	BOOST_CHECK( (*p.first >= -0.5) and (*p.second < 0.5) );

}
