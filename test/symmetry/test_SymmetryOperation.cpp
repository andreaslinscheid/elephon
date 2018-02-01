/*	This file test_SymmetryOperation.cpp is part of elephon.
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
 *  Created on: Jan 16, 2018
 *      Author: A. Linscheid
 */

#include <boost/test/unit_test.hpp>
#include "symmetry/SymmetryOperation.h"
#include "AtomicSite/ASSymmetry.h"
#include "fixtures/scenarios.h"
#include <cstdlib>
#include <ctime>
#include <cmath>

BOOST_AUTO_TEST_SUITE( SymmetryOperation )

BOOST_AUTO_TEST_CASE( test_creation )
{
	using namespace elephon;

	std::vector<int> ref_sym_pg{	0,1,0,
									0,0,1,
									1,0,0	};

	std::vector<double> ref_sym_ft{	0.0, 0.0, 0.0 };

	std::vector<double> ref_sym_pgc{	0,1,0,
										0,0,1,
										1,0,0	};

	std::vector<double> ref_sym_ftc{	0.0, 0.0, 0.0 };

	symmetry::SymmetryOperation sop(ref_sym_pg.begin(), ref_sym_pg.end(),
									ref_sym_ft.begin(), ref_sym_ft.end(),
									ref_sym_pgc.begin(), ref_sym_pgc.end(),
									ref_sym_ftc.begin(), ref_sym_ftc.end(),
									std::make_shared<AtomicSite::ASSymmetry::RadSym>());

	for (int i = 0 ; i < 3 ; ++i)
	{
		BOOST_CHECK_SMALL(ref_sym_ft[i]-sop.get_lat_frac_trans(i), 1e-5);
		BOOST_CHECK_SMALL(ref_sym_ftc[i]-sop.get_carth_frac_trans(i), 1e-5);
		for (int j = 0 ; j < 3 ; ++j)
		{
			BOOST_CHECK_EQUAL(ref_sym_pg[i*3+j], sop.get_lat_rot_matrix(i,j));
			BOOST_CHECK_SMALL(ref_sym_pgc[i*3+j] - sop.get_carth_rot_matrix(i,j), 1e-5);
		}
	}
}

BOOST_AUTO_TEST_CASE( test_apply_rotate )
{
	using namespace elephon;

	auto res = elephon::test::fixtures::scenarios::load_Al_fcc_primitive_vasp_sc4x4x4();
	auto sym = res->get_primitive_unitcell_obj()->get_symmetry();
	BOOST_REQUIRE( sym.get_num_symmetries() == 48 );

	// The symmetry operation 45 is the following matrix which exchanges x and y
	// 0, 1, 0
	// 1, 0, 0,
	// 0, 0, 1
	symmetry::SymmetryOperation sop = sym.get_sym_op(/*isym =*/45);

	std::vector<double> vec(3);
	srand(std::time(NULL));
	for (auto &v : vec)
	{
		v = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
		v -= std::floor(v+0.5); // map back 1. unit cell
	}
	auto vec_copy = vec;
	sop.apply(vec);

	// confirm the operation
	BOOST_CHECK_SMALL(vec[1] - vec_copy[0], 1e-5);
	BOOST_CHECK_SMALL(vec[0] - vec_copy[1], 1e-5);
	BOOST_CHECK_SMALL(vec[2] - vec_copy[2], 1e-5);
}

BOOST_AUTO_TEST_SUITE_END()
