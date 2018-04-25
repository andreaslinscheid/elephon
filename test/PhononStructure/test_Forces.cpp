/*	This file test_Forces.cpp is part of elephon.
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
 *  Created on: Feb 25, 2018
 *      Author: A. Linscheid
 */
#include <boost/test/unit_test.hpp>
#include "PhononStructure/Forces.h"
#include "fixtures/scenarios.h"
#include "LatticeStructure/PrimitiveToSupercellConnection.h"
#include "LatticeStructure/AtomDisplacementCollection.h"

BOOST_AUTO_TEST_SUITE( Forces )

BOOST_AUTO_TEST_CASE( Al_4x4x4_read_forces )
{
	auto resource = elephon::test::fixtures::scenarios::load_Al_fcc_primitive_vasp_sc4x4x4();
	auto displ = resource->get_displmts_collection_obj();
	elephon::PhononStructure::Forces forces;
	forces.initialize(displ , resource->get_electronic_structure_interface());

	BOOST_CHECK_EQUAL(forces.get_num_total_irred_displacements(), displ->get_tota_num_irred_displacements());

	const int displAtomIndex = 0;
	auto f = forces.get_forces_for_atom(displAtomIndex,0);
	BOOST_CHECK_EQUAL(f.shape()[0], resource->get_supercell_obj()->get_atoms_list().size());
	BOOST_CHECK_EQUAL(f.shape()[1], 3);

	elephon::test::fixtures::DataLoader dl;
	auto force_ref = dl.get_reference_force_data_vasp_Al_sc4x4x4();
	BOOST_REQUIRE_EQUAL(force_ref.size(), f.size());
	for (int iatomSC = 0 ; iatomSC < f.shape()[0]; ++iatomSC )
		for (int ix = 0 ;ix < 3; ++ix)
			BOOST_CHECK_CLOSE(f[iatomSC][ix], force_ref[iatomSC*3+ix], 1e-8);

	// Check the force on a given atom an that it transforms correctly.
	elephon::Auxillary::Multi_array<double,3> expandedForces;;
	forces.site_symmetry_expand_data(
			resource->get_primitive_unitcell_obj()->get_site_symmetry(displAtomIndex),
			displAtomIndex,
			resource->get_primitive_supercell_connect_obj()->primitive_to_supercell_atom_index(displAtomIndex),
			resource->get_supercell_obj()->get_atoms_list(),
			expandedForces);

	// force on the atom number 63 with coordinates (-0.25 -0.25 -0.25) is F = (0.00009213, 0.00000246, -0.00009213)
	// Symmetry operation number 9 is ((-1 -1 -1), ( 0  1  0), ( 0  0  1)) which takes the atom to
	// the site (-0.25 -0.25 -0.25), i.e. the same site.
	const int checksite = 63;
	const int checkSOP = 9;
	std::array<double, 3> f_check = {f[checksite][0], f[checksite][1], f[checksite][2]};
	BOOST_CHECK_CLOSE(f_check[0], 0.00009213, 1e-8);
	BOOST_CHECK_CLOSE(f_check[1], 0.00000246, 1e-8);
	BOOST_CHECK_CLOSE(f_check[2],-0.00009213, 1e-8);
	auto sop_9 = resource->get_primitive_unitcell_obj()->get_site_symmetry(displAtomIndex).get_sym_op(checkSOP);
	resource->get_primitive_unitcell_obj()->get_site_symmetry(displAtomIndex).rotate_cartesian(checkSOP, f_check.begin(), f_check.end());
	BOOST_CHECK_CLOSE(f_check[0], expandedForces[9][checksite][0], 1e-8);
	BOOST_CHECK_CLOSE(f_check[1], expandedForces[9][checksite][1], 1e-8);
	BOOST_CHECK_CLOSE(f_check[2], expandedForces[9][checksite][2], 1e-8);
}

BOOST_AUTO_TEST_SUITE_END()
