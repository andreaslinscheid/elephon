/*	This file test_UnitCell.cpp is part of elephon.
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
 *  Created on: May 15, 2017
 *      Author: A. Linscheid
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Input_test
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "LatticeStructure/UnitCell.h"
#include "IOMethods/VASPInterface.h"
#include "fixtures/MockStartup.h"
#include "fixtures/DataLoader.h"

elephon::LatticeStructure::UnitCell
load_unit_cell_Al_vasp_conventional()
{
	test::fixtures::MockStartup ms;
	auto testd = ms.get_data_for_testing_dir() / "Al" / "vasp" / "conventional";
	std::string input = std::string()+
			"root_dir = "+testd.string()+"\n";
	elephon::IOMethods::InputOptions opts;
	ms.simulate_elephon_input(
			(testd / "infile").string(),
			input,
			opts);

	auto loader = std::make_shared<elephon::IOMethods::VASPInterface>(opts);
	elephon::LatticeStructure::UnitCell uc;
	loader->read_unit_cell(testd.string(), opts.get_gPrec(), uc);
	return uc;
}

BOOST_AUTO_TEST_CASE( Generate_Al_fcc_primitive_displacements )
{
	test::fixtures::MockStartup ms;
	auto testd = ms.get_data_for_testing_dir() / "Al" / "vasp" / "fcc_primitive";
	std::string input = std::string()+
			"root_dir = "+testd.string()+"\n";
	elephon::IOMethods::InputOptions opts;
	ms.simulate_elephon_input(
			(testd / "infile").string(),
			input,
			opts);

	auto loader = std::make_shared<elephon::IOMethods::VASPInterface>(opts);
	elephon::LatticeStructure::UnitCell uc;
	loader->read_unit_cell(testd.string(), 1e-6, uc);

	std::vector<elephon::LatticeStructure::AtomDisplacement> irreducibleDisplacements;
	std::vector<int> mapRedToIrred, mapSymRedToIrred;
	std::vector<std::vector<int>> mapIrredToRed, mapSymIrredToRed;
	uc.generate_displacements( 0.01,
			/*bool symmetricDisplacement = */ true,
			irreducibleDisplacements);

	BOOST_REQUIRE_EQUAL(irreducibleDisplacements.size(),1);
}

BOOST_AUTO_TEST_CASE( Build_Al_supercell )
{
	auto uc = load_unit_cell_Al_vasp_conventional();

	elephon::LatticeStructure::UnitCell trivialSupercell = uc.build_supercell( 1, 1, 1);

	BOOST_REQUIRE( trivialSupercell.get_atoms_list().size() == 4 );

	elephon::LatticeStructure::UnitCell supercell = uc.build_supercell( 2, 1, 1);

	BOOST_REQUIRE( supercell.get_atoms_list().size() == 8 );
}

BOOST_AUTO_TEST_CASE( Generate_Al_displacements )
{
	auto uc = load_unit_cell_Al_vasp_conventional();

	std::vector<elephon::LatticeStructure::AtomDisplacement> irreducibleDisplacements;
	std::vector<int> mapRedToIrred, mapSymRedToIrred;
	std::vector<std::vector<int>> mapIrredToRed, mapSymIrredToRed;
	uc.generate_displacements( 0.01,
			/*bool symmetricDisplacement = */ false,
			irreducibleDisplacements);

	//how to test this?
}


BOOST_AUTO_TEST_CASE( Load_Al_vasp_fcc_primitve )
{
	using namespace elephon;
	test::fixtures::MockStartup ms;
	auto rootDir = ms.get_data_for_testing_dir() / "Al" / "vasp" / "fcc_primitive" / "phonon_run";

	std::string content = std::string("root_dir=")+rootDir.string()+"\n";
	test::fixtures::DataLoader dl;
	auto loader = dl.create_vasp_loader( content );

	LatticeStructure::UnitCell uc;
	loader->read_unit_cell(rootDir.string(), loader->get_optns().get_gPrec(),uc);

	BOOST_REQUIRE_EQUAL(uc.get_site_symmetry( 0 ).get_num_symmetries(), 48);
}
