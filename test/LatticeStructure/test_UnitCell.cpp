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
#include "IOMethods/ReadVASPSymmetries.h"
#include "LatticeStructure/Symmetry.h"
#include "LatticeStructure/LatticeModule.h"
#include "IOMethods/ReadVASPPoscar.h"

BOOST_AUTO_TEST_CASE( Build_Al_supercell )
{
	boost::filesystem::path p(__FILE__);
	boost::filesystem::path dir = p.parent_path();
	std::string Al_test_poscar = std::string(dir.c_str())+"/../IOMethods/POSCAR_Al_test.dat";

	elephon::IOMethods::ReadVASPPoscar filerreader;
	filerreader.read_file( Al_test_poscar );

	elephon::LatticeStructure::LatticeModule lattice;
	lattice.initialize(filerreader.get_lattice_matrix());

	elephon::IOMethods::ReadVASPSymmetries symreader;
	symreader.read_file( std::string(dir.c_str())+"/../IOMethods/OUTCAR_Al_test.dat" );
	elephon::LatticeStructure::Symmetry sym;
	sym.initialize( 1e-6 , symreader.get_symmetries(), symreader.get_fractionTranslations());

	elephon::LatticeStructure::UnitCell uc;
	uc.initialize( filerreader.get_atoms_list(), lattice, sym);

	elephon::LatticeStructure::UnitCell trivialSupercell = uc.build_supercell( 1, 1, 1);

	BOOST_REQUIRE( trivialSupercell.get_atoms_list().size() == 4 );

	elephon::LatticeStructure::UnitCell supercell = uc.build_supercell( 2, 1, 1);

	BOOST_REQUIRE( supercell.get_atoms_list().size() == 8 );
}

BOOST_AUTO_TEST_CASE( Generate_Al_displacements )
{
	using namespace elephon;

	boost::filesystem::path p(__FILE__);
	boost::filesystem::path dir = p.parent_path();
	std::string Al_test_poscar = std::string(dir.c_str())+"/../IOMethods/POSCAR_Al_test.dat";

	IOMethods::ReadVASPPoscar filerreader;
	filerreader.read_file( Al_test_poscar );

	LatticeStructure::LatticeModule lattice;
	lattice.initialize(filerreader.get_lattice_matrix());

	IOMethods::ReadVASPSymmetries symreader;
	symreader.read_file( std::string(dir.c_str())+"/../IOMethods/OUTCAR_Al_test.dat" );
	LatticeStructure::Symmetry sym;
	sym.initialize( 1e-6 , symreader.get_symmetries(), symreader.get_fractionTranslations());

	LatticeStructure::UnitCell uc;
	uc.initialize( filerreader.get_atoms_list(), lattice, sym);

	std::vector<LatticeStructure::AtomDisplacement> reducibleDisplacements,irreducibleDisplacements;
	std::vector<int> mapRedToIrred, mapSymRedToIrred;
	std::vector<std::vector<int>> mapIrredToRed, mapSymIrredToRed;
	uc.generate_displacements( 0.01,
			/*bool symmetricDisplacement = */ false,
			reducibleDisplacements,irreducibleDisplacements,
			mapRedToIrred, mapSymRedToIrred,
			mapIrredToRed, mapSymIrredToRed);

	size_t numModes = 3*4;// #Atoms = 4
	BOOST_REQUIRE( reducibleDisplacements.size() == numModes);

	uc.generate_displacements( 0.01,
			/*bool symmetricDisplacement = */ true,
			reducibleDisplacements,irreducibleDisplacements,
			mapRedToIrred, mapSymRedToIrred,
			mapIrredToRed, mapSymIrredToRed);

	numModes = 2*3*4;// #Atoms = 4, +-xi
	BOOST_REQUIRE( reducibleDisplacements.size() == numModes);
}
