/*	This file test_RegularBareGrid.cpp is part of elephon.
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
 *  Created on: Sep 27, 2017
 *      Author: A. Linscheid
 */

#define BOOST_TEST_MODULE Input_test
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "LatticeStructure/LatticeModule.h"
#include "LatticeStructure/RegularBareGrid.h"
#include "IOMethods/ReadVASPPoscar.h"
#include "IOMethods/ReadVASPKpoints.h"
#include "fixtures/MockStartup.h"
#include <vector>

BOOST_AUTO_TEST_CASE( G2S2_vasp_bare_kpoint_grid )
{
	elephon::test::fixtures::MockStartup ms;
	auto path = ms.get_data_for_testing_dir() / "Ge2S2" ;

	std::vector<std::pair<std::string, double> > atoms;
	// mass Ge 72.63
	// mass S 32.065
	atoms.push_back(std::make_pair("Ge", 72.63) );
	atoms.push_back(std::make_pair("S", 32.065) );

	elephon::IOMethods::ReadVASPPoscar poscarReader;
	poscarReader.read_file( (path / "POSCAR").string(), atoms );

	elephon::LatticeStructure::LatticeModule lattice;
	lattice.initialize( poscarReader.get_lattice_matrix() );
	elephon::IOMethods::ReadVASPKpoints kpointsReader;
	kpointsReader.read_kpoints(  (path / "KPOINTS").string(), lattice );
	elephon::LatticeStructure::RegularBareGrid kgrid;
	kgrid.initialize(kpointsReader.get_grid_dim(), true, 1e-6, kpointsReader.get_grid_shift(), lattice);
}

