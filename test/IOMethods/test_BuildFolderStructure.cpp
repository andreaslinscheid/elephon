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
#include "fixtures/MockStartup.h"
#include <vector>

BOOST_AUTO_TEST_CASE( Build_Al_primitive_folderstructure_VASP )
{
	using namespace boost::filesystem;

	test::fixtures::MockStartup ms;
	auto testd = ms.get_data_for_testing_dir();
	elephon::IOMethods::BuildFolderStructure builder;

	path rootDir = testd / "Al" / "vasp" / "fcc_primitive";
	path targetDir = rootDir / "phonon_test_build";

	//here we create the test input file
	std::string content = std::string()+
			"scell=2 2 2\n"
			"root_dir="+rootDir.string()+"\n"
			"elphd="+targetDir.string()+"\n"
			"";
	std::string filename = (rootDir / "test_elephon_input.dat").string();

	remove_all( targetDir  );

	elephon::IOMethods::InputOptions options;
	ms.simulate_elephon_input( (rootDir / "test_folderstructure_input.dat").string(),
			content, options );

	//Here we use the VASP interface to generate the folder structure
	elephon::IOMethods::VASPInterface vi(options);
	elephon::LatticeStructure::Symmetry sym;
	std::vector<elephon::LatticeStructure::Atom> atoms;
	elephon::LatticeStructure::LatticeModule lattice;
	elephon::LatticeStructure::RegularSymmetricGrid kgrid;
	vi.read_cell_paramters(
			rootDir.string(),
			1e-6,
			kgrid,
			lattice,
			atoms,
			sym);

	elephon::LatticeStructure::UnitCell uc;
	uc.initialize( atoms, lattice, sym);

	builder.build( options, uc, vi );

	BOOST_REQUIRE( builder.check_is_build( targetDir.string() )  ==  true );

	BOOST_REQUIRE( exists(targetDir / "electrons" )  == true );
	BOOST_REQUIRE( exists(targetDir / "electrons" / "POSCAR" )  == true );
	BOOST_REQUIRE( exists(targetDir / "electrons" / "KPOINTS" )  == true );
	BOOST_REQUIRE( exists(targetDir / "electrons" / "INCAR" )  == true );
	BOOST_REQUIRE( exists(targetDir / "electrons" / "POTCAR" )  == true );

	//We have 1 irreducible displacements
	BOOST_REQUIRE( exists(targetDir / "displ_0" )  == true );
	BOOST_REQUIRE( exists(targetDir / "displ_0" / "POSCAR" )  == true );
	BOOST_REQUIRE( exists(targetDir / "displ_0" / "KPOINTS" )  == true );
	BOOST_REQUIRE( exists(targetDir / "displ_0" / "INCAR" )  == true );
	BOOST_REQUIRE( exists(targetDir / "displ_0" / "POTCAR" )  == true );
	BOOST_REQUIRE( exists(targetDir / "displ_1" )  == false );
}

BOOST_AUTO_TEST_CASE( Build_Al_folderstructure_VASP )
{
	using namespace boost::filesystem;
	elephon::IOMethods::BuildFolderStructure builder;
	test::fixtures::MockStartup ms;
	auto testd = ms.get_data_for_testing_dir();
	path rootDir = testd / "Al" / "vasp" / "conventional";
	path targetDir = rootDir / "test_dir_struct";

	//here we create the test input file
	std::string content = std::string()+
			"scell=1 1 1\n"
			"root_dir="+rootDir.string()+"\n"
			"elphd="+targetDir.string()+"\n"
			"";
	std::string filename = (rootDir / "test_elephon_input.dat").string();

	elephon::IOMethods::InputOptions options;
	ms.simulate_elephon_input( (rootDir / "test_folderstructure_input.dat").string(),
			content, options );

	//Here we use the VASP interface to generate the folder structure
	elephon::IOMethods::VASPInterface vi(options);
	elephon::LatticeStructure::Symmetry sym;
	std::vector<elephon::LatticeStructure::Atom> atoms;
	elephon::LatticeStructure::LatticeModule lattice;
	elephon::LatticeStructure::RegularSymmetricGrid kgrid;
	vi.read_cell_paramters(
			rootDir.string(),
			1e-6,
			kgrid,
			lattice,
			atoms,
			sym);

	elephon::LatticeStructure::UnitCell uc;
	uc.initialize( atoms, lattice, sym);

	remove_all(targetDir);

	BOOST_CHECK( builder.check_is_build( targetDir.string() )  == false );

	builder.build( options, uc, vi );

	BOOST_REQUIRE( builder.check_is_build( targetDir.string() )  ==  true );

	BOOST_REQUIRE( exists(targetDir / "electrons" )  == true );

	//We have 4 irreducible displacements
	BOOST_REQUIRE( exists(targetDir / "displ_0" )  == true );
	BOOST_REQUIRE( exists(targetDir / "displ_1" )  == true );
	BOOST_REQUIRE( exists(targetDir / "displ_2" )  == true );
	BOOST_REQUIRE( exists(targetDir / "displ_3" )  == true );
	BOOST_REQUIRE( exists(targetDir / "displ_4" )  == false );
}
