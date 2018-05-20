/*	This file test_EliashbergModule.cpp is part of elephon.
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
 *  Created on: May 18, 2018
 *      Author: A. Linscheid
 */

#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "EliashbergEquations/EliashbergModule.h"
#include "fixtures/MockStartup.h"
#include "fixtures/DataLoader.h"

BOOST_AUTO_TEST_SUITE( EliashbergModule )

BOOST_AUTO_TEST_CASE( test_Tc_zero_coulomb )
{
	elephon::test::fixtures::MockStartup ms;
	auto testd = ms.get_data_for_testing_dir() / "MgB2" / "ref_data";
	elephon::test::fixtures::DataLoader dl;
	auto res = dl.create_resource_handler(std::string("")+
			"f_a2F = "+(testd / "a2F_isotropic_THz.dat").string()+"\n"
			"EliT = findTc\n"
			"muStar = 0.0\n");
	elephon::EliashbergEquations::EliashbergModule eliMod(res->get_optns());
	eliMod.do_work();
}

BOOST_AUTO_TEST_SUITE_END()
