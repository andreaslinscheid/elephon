/*	This file test_ReadVASPElephonData.cpp is part of elephon.
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
 *  Created on: Nov 27, 2017
 *      Author: A. Linscheid
 */
#include <boost/test/unit_test.hpp>
#include "fixtures/scenarios.h"

BOOST_AUTO_TEST_SUITE( ReadVASPElephonData )

BOOST_AUTO_TEST_CASE( Read_Al_fcc_primitive_vasp_data )
{
	auto resHandler = elephon::test::fixtures::scenarios::load_Al_fcc_primitive_vasp_sc2x2x2();
//	auto atoms = resHandler->get_primitive_unitcell_obj()->get_atoms_list();
//
//	BOOST_REQUIRE(atoms.size() == 1);
//	auto pawModule = atoms[0].get_paw_data();
//	BOOST_REQUIRE(pawModule == false);
//
//	resHandler->load_PAW_data();
//
//	pawModule = resHandler->get_primitive_unitcell_obj()->get_atoms_list()[0].get_paw_data();
//	BOOST_REQUIRE(pawModule);
}

BOOST_AUTO_TEST_SUITE_END()
