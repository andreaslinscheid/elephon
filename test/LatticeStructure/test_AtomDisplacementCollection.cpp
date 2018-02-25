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

	BOOST_REQUIRE_EQUAL(displ->get_tota_num_irred_displacements(),1);
	BOOST_REQUIRE_EQUAL(displ->get_tota_num_irred_displacements(), displ->get_irreducible_displacements().size());
	BOOST_CHECK_EQUAL(displ->get_irreducible_displacements()[0].second.size(),1);

	// Atom no. 0 is at a high symmetry position (0,0,0). There are 48 symmetries, i.e. subtracting inversion, there are
	BOOST_CHECK_EQUAL(displ->get_red_displacements_for_atom(0).size(), 48/2 );
}
BOOST_AUTO_TEST_SUITE_END()
