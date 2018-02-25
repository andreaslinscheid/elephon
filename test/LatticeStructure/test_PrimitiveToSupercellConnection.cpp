/*	This file test_PrimitiveToSupercellConnection.cpp is part of elephon.
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
 *  Created on: Feb 05, 2017
 *      Author: A. Linscheid
 */
#include <boost/test/unit_test.hpp>
#include "LatticeStructure/PrimitiveToSupercellConnection.h"
#include "fixtures/scenarios.h"

BOOST_AUTO_TEST_SUITE( PrimitiveToSupercellConnection )

BOOST_AUTO_TEST_CASE( EqualScaleEveryDirection )
{
	auto resource = elephon::test::fixtures::scenarios::load_Al_fcc_primitive_vasp_sc4x4x4();
	auto primitiveCell = resource->get_primitive_unitcell_obj();
	auto superCell = resource->get_supercell_obj();

	elephon::LatticeStructure::PrimitiveToSupercellConnection primSCCon;
	primSCCon.initialize( primitiveCell, superCell );

	const int scX = 4;
	const int scY = 4;
	const int scZ = 4;
	BOOST_CHECK( primSCCon.get_supercell_volume_factor() == scX*scY*scZ );

	auto AlAtom = primitiveCell->get_atoms_list().front();
	BOOST_REQUIRE(AlAtom.get_kind().compare("Al")==0);
	auto index = primSCCon.find_atom(AlAtom, true);
	// The following checks AlAtom == primitiveCell->get_atoms_list()[index]
	BOOST_CHECK(AlAtom == primitiveCell->get_atoms_list()[index]);

	// the atom in the primitve cell is at 0. We check if we get the atom at
	srand(time(nullptr));
	const int Rx = rand() % scX;
	const int Ry = rand() % scY;
	const int Rz = rand() % scZ;
	BOOST_MESSAGE("Got vector "+std::to_string(Rx)+", "+std::to_string(Ry)+", "+std::to_string(Rz));

	std::vector<double> checkVector{ 	static_cast<double>(Rx)/static_cast<double>(scX),
										static_cast<double>(Ry)/static_cast<double>(scY),
										static_cast<double>(Rz)/static_cast<double>(scZ) };
	AlAtom.set_position(checkVector);
	index = primSCCon.find_atom(AlAtom, false);
	BOOST_CHECK(AlAtom == superCell->get_atoms_list()[index]);

	auto oldVector = checkVector;
	primSCCon.supercell_to_primitive_coordinates(checkVector);
	BOOST_CHECK_CLOSE(checkVector[0], static_cast<double>(Rx),1e-6);
	BOOST_CHECK_CLOSE(checkVector[1], static_cast<double>(Ry),1e-6);
	BOOST_CHECK_CLOSE(checkVector[2], static_cast<double>(Rz),1e-6);

	primSCCon.primitive_to_supercell_coordinates(checkVector);
	BOOST_CHECK_CLOSE(checkVector[0], oldVector[0], 1e-6);
	BOOST_CHECK_CLOSE(checkVector[1], oldVector[1], 1e-6);
	BOOST_CHECK_CLOSE(checkVector[2], oldVector[2], 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()
