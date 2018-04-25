/*	This file test_AtomDisplacementCollection.cpp is part of elephon.
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
 *  Created on: Feb 20, 2018
 *      Author: A. Linscheid
 */
#include <boost/test/unit_test.hpp>
#include "LatticeStructure/AtomDisplacementCollection.h"
#include "fixtures/MockStartup.h"
#include "fixtures/DataLoader.h"

BOOST_AUTO_TEST_SUITE( AtomDisplacementCollection )

BOOST_AUTO_TEST_CASE( Generate_Al_fcc_primitive_displacements )
{
	elephon::test::fixtures::MockStartup ms;
	auto testd = ms.get_data_for_testing_dir() / "Al" / "vasp" / "fcc_primitive";
	std::string input = std::string()+
			"root_dir = "+testd.string()+"\n";
	elephon::test::fixtures::DataLoader dl;
	auto resource =  dl.create_resource_handler(input);
	// this actually generates the displacements
	auto displ = resource->get_displmts_collection_obj();
	auto primitiveCell = resource->get_primitive_unitcell_obj();

	BOOST_REQUIRE_EQUAL(displ->get_tota_num_irred_displacements(),1);
	BOOST_REQUIRE_EQUAL(displ->get_tota_num_irred_displacements(), displ->get_irreducible_displacements().size());
	BOOST_CHECK_EQUAL(displ->get_irreducible_displacements()[0].second.size(),1);

	// Atom no. 0 is at a high symmetry position (0,0,0). There are 48 symmetries there are
	const int AtomIndex = 0;
	auto reducibleDisplacements = displ->get_red_displacements_for_atom(AtomIndex);
	BOOST_CHECK_EQUAL(reducibleDisplacements.size(), 48 );

	BOOST_CHECK(displ->check_atom_is_irreducible(0)); // Only one atom which is irreducible
	BOOST_CHECK(displ->get_num_atoms_primitive_cell() == 1);

	auto atomAndRelIrred = displ->get_total_irred_index_to_atom_and_rel_irred_index(0);
	BOOST_CHECK_EQUAL(atomAndRelIrred.first,0); // First Atom
	BOOST_CHECK_EQUAL(atomAndRelIrred.second,0); //First irreducible displacement

	auto redDisplSymRel = displ->get_symmetry_relation_red_displacements_for_atom(AtomIndex);
	BOOST_REQUIRE_EQUAL(redDisplSymRel.size(), reducibleDisplacements.size());
	auto siteSymmetry = primitiveCell->get_site_symmetry(AtomIndex);
	auto irrd = displ->get_irreducible_displacements();
	BOOST_REQUIRE_EQUAL(irrd.size(), 1);
	BOOST_REQUIRE_EQUAL(irrd[0].second.size(), 1);
	elephon::Auxillary::Multi_array<double,2> displacements(boost::extents[48][3]);
	for (int redDisplIndex = 0; redDisplIndex < redDisplSymRel.size(); ++redDisplIndex )
	{
		BOOST_CHECK_EQUAL(redDisplSymRel[redDisplIndex].first, 0); // Only one irreducible displacement
		const int symIDSiteSymIrredToRed = redDisplSymRel[redDisplIndex].second;
		auto sop = siteSymmetry.get_sym_op(symIDSiteSymIrredToRed);
		auto thisIrredDispl = irrd[0].second[0];
		thisIrredDispl.transform(sop);
		// check if the transformed irreducible displacement equals the reducible one, as required
		BOOST_CHECK(thisIrredDispl == reducibleDisplacements[redDisplIndex]);
		for (int i = 0 ; i < 3; ++i)
			displacements[redDisplIndex][i] = reducibleDisplacements[redDisplIndex].get_magnitude()*
					reducibleDisplacements[redDisplIndex].get_direction()[i];
	}

	// here is the pseudo inverse of the displacement matrix, the function is supposed to return
	elephon::Algorithms::LinearAlgebraInterface linalg;
	elephon::Auxillary::Multi_array<double,2> pseudoInvU_ref(boost::extents[48][3]), pseudoInvU;
	linalg.pseudo_inverse(displacements, 48, 3, pseudoInvU_ref);
	pseudoInvU_ref.reshape(boost::array<int,2>({{3,48}}));
	displ->generate_pseudo_inverse_reducible_displacements_for_atom(AtomIndex, pseudoInvU);
	BOOST_REQUIRE_EQUAL(pseudoInvU.shape()[0], pseudoInvU_ref.shape()[0]);
	BOOST_REQUIRE_EQUAL(pseudoInvU.shape()[1], pseudoInvU_ref.shape()[1]);
	double sumDiff = 0.0;
	for (int i = 0 ; i < pseudoInvU.shape()[0]; ++i)
		for (int j = 0 ; j < pseudoInvU.shape()[1]; ++j)
			sumDiff += std::abs(pseudoInvU_ref[i][j]-pseudoInvU[i][j]);
	const double tolerance = 1e-8;
	BOOST_CHECK(sumDiff<tolerance);
}
BOOST_AUTO_TEST_SUITE_END()
