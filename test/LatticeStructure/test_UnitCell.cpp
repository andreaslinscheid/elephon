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

	elephon::LatticeStructure::Symmetry sym;

	elephon::LatticeStructure::UnitCell uc;
	uc.initialize( filerreader.get_atoms_list(), lattice, sym);

	elephon::LatticeStructure::UnitCell supercell = uc.build_supercell( 2, 1, 1);

	BOOST_REQUIRE( supercell.get_atoms_list().size() == 8 );
}
