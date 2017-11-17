/*	This file test_ReadVASPPotcar.cpp is part of elephon.
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
 *  Created on: May 31, 2017
 *      Author: A. Linscheid
 */

#define BOOST_TEST_MODULE Input_test
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include <IOMethods/ReadVASPLocpot.h>
#include "fixtures/MockStartup.h"
#include <assert.h>
#include <vector>

BOOST_AUTO_TEST_CASE( Read_potcar )
{
	elephon::IOMethods::ReadVASPLocpot filerreader;
	elephon::test::fixtures::MockStartup ms;
	auto testd = ms.get_data_for_testing_dir() / "Al" / "vasp" / "conventional" / "supercell";
	std::vector<double> potential;
	std::vector<int> dims;
	filerreader.read_scf_potential(
			(testd / "LOCPOT").string(),
			dims,
			potential);

	std::vector<int> fftDims = {60,60,30};

	BOOST_REQUIRE( fftDims == dims );

	BOOST_REQUIRE( int(potential.size()) == fftDims[0]*fftDims[1]*fftDims[2] );

	double firstVal = -.12739458907E+03;
	double butLastVal = -.10182458035E+03;

	BOOST_REQUIRE( std::abs(potential[0] - firstVal) < 1e-8 );
	double testVal = *(++potential.rbegin());
	BOOST_REQUIRE( std::abs( testVal - butLastVal) < 1e-8 );
}
