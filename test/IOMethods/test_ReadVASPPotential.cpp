/*	This file test_ReadVASPPotential.cpp is part of elephon.
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
 *  Created on: Apr 26, 2018
 *      Author: A. Linscheid
 */

#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "IOMethods/ReadVASPPotential.h"
#include "fixtures/MockStartup.h"
#include <assert.h>
#include <vector>
#include <array>

BOOST_AUTO_TEST_SUITE( ReadVASPPotential )

BOOST_AUTO_TEST_CASE( Read_Al_LOCPOT_AE )
{
	elephon::IOMethods::ReadVASPPotential filerreader;
	elephon::test::fixtures::MockStartup ms;
	filerreader.set_filepath((ms.get_data_for_testing_dir() / "Al" / "vasp" / "fcc_primitive" / "sc_4x4x4"
			/ filerreader.get_default_filename()).string());

	std::array<int,3> regularGridDims;
	std::vector<double> regularData, radiusPerAtom, coreChargeZ;
	std::vector<int> angularL;
	std::vector<std::vector<double>> radialGrid, coreFrozenElectronicCharge;
	std::vector<std::vector<std::complex<double>>> radialPotentialDatal;
	filerreader.read_potential_file(regularGridDims, regularData, angularL,
			radialGrid, radiusPerAtom, radialPotentialDatal, coreChargeZ, coreFrozenElectronicCharge);

	BOOST_REQUIRE_EQUAL( radialGrid.size(), 1u ); // one Atom in this unit cell
	BOOST_REQUIRE_EQUAL( radialGrid.size(), angularL.size() );
	BOOST_REQUIRE_EQUAL( radialGrid.size(), radiusPerAtom.size() );
	BOOST_REQUIRE_EQUAL( radialGrid.size(), radialPotentialDatal.size() );

	std::array<int,3> fftDims = {32,32,32};

	BOOST_REQUIRE( regularGridDims == fftDims );

	BOOST_REQUIRE_EQUAL( std::int64_t(regularData.size()), fftDims[0]*fftDims[1]*fftDims[2] );

	double firstVal = -.12760769678E+03;
	double butLastVal = -.10585117804E+03;

	BOOST_CHECK_SMALL( std::abs(regularData[0] - firstVal), 1e-5 );
	double testVal = *(++regularData.rbegin());
	BOOST_CHECK_SMALL( std::abs( testVal - butLastVal), 1e-5 );
}

BOOST_AUTO_TEST_SUITE_END()
