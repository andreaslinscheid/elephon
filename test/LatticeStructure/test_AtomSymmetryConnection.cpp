/*	This file test_AtomSymmetryConnection.cpp is part of elephon.
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
 *  Created on: Feb 27, 2018
 *      Author: A. Linscheid
 */
#include <boost/test/unit_test.hpp>
#include "fixtures/scenarios.h"
#include "LatticeStructure/AtomSymmetryConnection.h"

BOOST_AUTO_TEST_SUITE( AtomSymmetryConnection )

BOOST_AUTO_TEST_CASE( Al_trivial )
{
	auto res = elephon::test::fixtures::scenarios::load_Al_fcc_primitive_vasp_sc2x2x2();
	auto atomSymmetry = res->get_primitive_unitcell_obj()->get_atom_symmetry();

	const int idIndex = res->get_primitive_unitcell_obj()->get_symmetry().get_identity_index();
	const int AlAtomIndex = 0;
	BOOST_CHECK(atomSymmetry->check_atom_is_irreducible(AlAtomIndex));
	auto star = atomSymmetry->get_star_atom_indices(AlAtomIndex);
	BOOST_CHECK(star.size() == 1);
	BOOST_CHECK(star[0].first == AlAtomIndex);
	BOOST_CHECK(star[0].second == AlAtomIndex);
	BOOST_CHECK(atomSymmetry->atom_rot_map(idIndex, AlAtomIndex) == AlAtomIndex);
	BOOST_CHECK(atomSymmetry->get_list_irreducible_atoms().size() == 1);
	BOOST_CHECK(atomSymmetry->get_list_irreducible_atoms()[0] == AlAtomIndex);
}

BOOST_AUTO_TEST_SUITE_END()
