/*	This file test_BuildFolderStructure.cpp is part of elephon.
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
 *  Created on: May 17, 2017
 *      Author: A. Linscheid
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Input_test
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "IOMethods/BuildFolderStructure.h"
#include "IOMethods/VASPInterface.h"
#include "IOMethods/Input.h"
#include "IOMethods/InputOptions.h"
#include <vector>

BOOST_AUTO_TEST_CASE( Build_Al_folderstructure_VASP )
{
	using namespace boost::filesystem;
	elephon::IOMethods::BuildFolderStructure builder;

	path dir = path(__FILE__).parent_path();
	path test_elph_dir = dir / "test_dir_struct";

	//here we create the test input file
	path test_input_file = dir / "test_folderstructure_input.dat";
	std::string content = std::string()+
			"scell=2 2 1\n"
			"root_dir="+(dir/"Al_test").string()+"\n"
			"elphd="+test_elph_dir.string()+"\n"
			"";
	std::ofstream file( test_input_file.c_str() );
	file << content;
	file.close();

	//here we reed the input file via elephons input mechanism
	char * prog = strdup("program name");
	char * arg = strdup(test_input_file.c_str());
	char *argv[] = {prog, arg, NULL};
	int argc = sizeof(argv) / sizeof(char*) - 1;
	elephon::IOMethods::Input input(argc,argv);
	elephon::IOMethods::InputOptions options = input.get_opts();
	delete [] prog;
	delete [] arg;

	//Here we use the VASP interface to generate the folder structure
	elephon::IOMethods::VASPInterface vi(options);
	elephon::LatticeStructure::Symmetry sym;
	std::vector<elephon::LatticeStructure::Atom> atoms;
	elephon::LatticeStructure::LatticeModule lattice;
	elephon::LatticeStructure::RegularGrid kgrid;
	vi.read_cell_paramters(
			(dir / "Al_test").string(),
			1e-6,
			kgrid,
			lattice,
			atoms,
			sym);

	elephon::LatticeStructure::UnitCell uc;
	uc.initialize( atoms, lattice, sym);

	remove_all(test_elph_dir);

	BOOST_CHECK( builder.check_is_build( test_elph_dir.string() )  == false );

	builder.build( options, uc, vi );

	BOOST_REQUIRE( builder.check_is_build( test_elph_dir.string() )  ==  true );

	BOOST_REQUIRE( exists(test_elph_dir / "electrons" )  == true );

	BOOST_REQUIRE( exists(test_elph_dir / "scell" )  == true );

	//We have 18 irreducible displacements
	BOOST_REQUIRE( exists(test_elph_dir / "displ_0" )  == true );
	BOOST_REQUIRE( exists(test_elph_dir / "displ_1" )  == true );
	BOOST_REQUIRE( exists(test_elph_dir / "displ_2" )  == true );
	BOOST_REQUIRE( exists(test_elph_dir / "displ_17" )  == true );
	BOOST_REQUIRE( exists(test_elph_dir / "displ_18" )  == false );
}
