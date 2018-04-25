/*	This file test_helperfunctions.cpp is part of elephon.
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
 *  Created on: Feb 12, 2018
 *      Author: A. Linscheid
 */
#include <boost/test/unit_test.hpp>
#include "Algorithms/helperfunctions.hpp"
#include "LatticeStructure/RegularBareGrid.h"
#include "Algorithms/GridRotationMap.h"
#include "fixtures/DataLoader.h"
#include <vector>

BOOST_AUTO_TEST_SUITE( helperfunctions )

template<typename T>
void
test_transform_vector_field(
		T xVal,
		elephon::LatticeStructure::RegularBareGrid const & grid,
		elephon::symmetry::SymmetryOperation const & sop,
		std::vector<int> const & map )
{
	std::vector<T> field(grid.get_num_points()*3, T(0));
	for (int i = 0 ; i < grid.get_num_points(); ++i)
		field[i*3] = (xVal*i)/grid.get_grid_dim()[0];

	std::vector<T> matrix(9);
	for (int i = 0 ; i < 3; ++i)
		for (int j = 0 ; j < 3; ++j)
			matrix[i*3+j] = sop.get_carth_rot_matrix(i,j);
	elephon::Algorithms::helperfunctions::transform_vector_field_cart(field.begin(), field.end(), matrix, map);

	double sumDiff = 0.0;
	for (int i = 0 ; i < grid.get_num_points(); ++i)
		sumDiff += (field[map[i]*3+1] + (xVal*i)/grid.get_grid_dim()[0]);
	BOOST_CHECK_SMALL(sumDiff, 1e-5);
}

BOOST_AUTO_TEST_CASE( transform_vector_field )
{
	// define a vector field with arrows pointing in the x direction with a magnitude proportional to x,
	// then apply a rotation by 90 Deg where the arrows are now  supposed to point in the y direction, magnitude prop to y.
	elephon::LatticeStructure::RegularBareGrid grid({10,10,10});
	elephon::test::fixtures::DataLoader dl;
	elephon::LatticeStructure::Symmetry sym = dl.create_partial_sym(); // sym no. 2 is the 90deg rotation in question
	std::vector<std::vector<int>> rotmap;
	elephon::Algorithms::compute_grid_rotation_map_no_shift(grid, sym, rotmap);

	// by explicitely checking that r=(1,0,0) ~ ir == 1, the first non-zero index in x is mapped
	// correctly, we convince ourself that this the correct mapping.
	const int indexRotatedFirstX1 = grid.get_xyz_to_reducible({0,1,0});
	BOOST_REQUIRE_EQUAL(indexRotatedFirstX1, rotmap[2][1]);

	test_transform_vector_field<float>(4.0f, grid, sym.get_sym_op(2), rotmap[2]);
}

BOOST_AUTO_TEST_SUITE_END()
