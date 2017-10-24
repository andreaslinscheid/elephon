/*	This file test_NestingFunction.cpp is part of elephon.
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
 *  Created on: Oct 22, 2017
 *      Author: A. Linscheid
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE NestingFunction
#include <boost/test/unit_test.hpp>
#include "ElectronicStructure/NestingFunction.h"
#include "IOMethods/ResourceHandler.h"
#include "LatticeStructure/TetrahedraGrid.h"
#include "ElectronicStructure/TetrahedraIsosurface.h"
#include "fixtures/MockStartup.h"
#include "fixtures/DataLoader.h"
#include <vector>

BOOST_AUTO_TEST_CASE( Al_fcc_vasp_nesting_function )
{
	elephon::test::fixtures::MockStartup ms;
	auto rootDir = ms.get_data_for_testing_dir() / "Al" / "vasp" / "fcc_primitive" ;
	elephon::test::fixtures::DataLoader dl;
	auto resourceHandler = dl.create_resource_handler( std::string()+
			"root_dir="+rootDir.string()+"\n");

	auto bands = resourceHandler->get_electronic_bands_obj();

	auto tetra = std::make_shared<elephon::LatticeStructure::TetrahedraGrid>();
	tetra->initialize( std::make_shared<elephon::LatticeStructure::RegularSymmetricGrid>(bands->get_grid()) );

	auto tetraIso = std::make_shared<elephon::ElectronicStructure::TetrahedraIsosurface>();
	tetraIso->initialize(tetra, bands, {0.0});
	auto grid = tetraIso->get_tetra_grid()->get_grid();

	elephon::ElectronicStructure::NestingFunction nest;
	nest.compute_nesting_function( tetraIso );

	auto isoNest = nest.get_nesting(0);

	std::vector<double> dosEf;
	bands->compute_DOS_tetra( tetra, {0.0}, dosEf);

	// integrate the nesting function
	double integralNest = 0;
	double dV = grid->get_lattice().get_reci_volume()
					/ grid->get_np_red();
	for ( int iq = 0 ; iq < grid->get_np_irred() ; ++iq)
	{
		int multi = grid->get_maps_sym_irred_to_reducible()[iq].size();
		for ( int ib = 0 ; ib < isoNest->get_nData_gpt() ; ++ib)
			integralNest += isoNest->read(iq, ib)*multi*dV;
	}

	 // note that this is not a very accurate calculation due to the
	 // problematic delta function approximation of the q grid
	std::cout << " Comparison of integrated nesting function and DOS squared: "
			<< integralNest << "(integral) vs " << dosEf[0]*dosEf[0] << "(DOS^2)"  << std::endl;
	BOOST_CHECK_CLOSE(integralNest, dosEf[0]*dosEf[0], 20 );
}
