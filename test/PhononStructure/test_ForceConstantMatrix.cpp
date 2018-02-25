/*	This file test_ForceConstantMatrix.cpp is part of elephon.
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
 *  Created on: Jun 1, 2017
 *      Author: A. Linscheid
 */
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "PhononStructure/ForceConstantMatrix.h"
#include "PhononStructure/Forces.h"
#include "LatticeStructure/AtomDisplacementCollection.h"
#include "IOMethods/Input.h"
#include "IOMethods/VASPInterface.h"
#include "fixtures/MockStartup.h"
#include "fixtures/DataLoader.h"
#include "fixtures/scenarios.h"
#include <assert.h>
#include <vector>
#include <fstream>
#include <string>
#include <memory>

BOOST_AUTO_TEST_SUITE( ForceConstantMatrix )

BOOST_AUTO_TEST_CASE( ForceConstantMatrix_symmetrization )
{
	std::vector<double> q = {0.5, 0.0, 0.0};
	std::complex<double> ii(0,1);
	std::complex<double> one(1,0);
	elephon::Auxillary::alignedvector::ZV dataV{	one, ii, one,
													one, ii, one,
													one, ii, one };
	elephon::Auxillary::Multi_array<std::complex<double>,3> data(boost::extents[1][3][3]);
	std::copy_n(dataV.begin(), dataV.size(), data.data());
	auto resHndl = elephon::test::fixtures::scenarios::load_Al_fcc_primitive_vasp_sc2x2x2();
	auto uc = resHndl->get_primitive_unitcell_obj();
	auto phi = resHndl->get_forceConstant_obj();
	phi->symmetrize_q(q, data, uc);
}

BOOST_AUTO_TEST_CASE( build_Al_primitive )
{
	elephon::test::fixtures::MockStartup ms;
	auto rootDir = ms.get_data_for_testing_dir() / "Al" / "vasp" / "fcc_primitive" ;
	auto phononDir = rootDir / "el_ph";
	std::string content = std::string()+
			"scell=2 2 2\n"
			"root_dir="+rootDir.string()+"\n"
			"elphd="+phononDir.string()+"\n"
			"";
	elephon::test::fixtures::DataLoader dl;
	auto res = dl.create_resource_handler(content);
	auto primitiveCell = res->get_primitive_unitcell_obj();
	auto supercell = res->get_supercell_obj();
	auto displColl = res->get_displmts_collection_obj();
	auto primToSCConn = res->get_primitive_supercell_connect_obj();

	//Here, we read the forces from the vasp output
	auto forces = std::make_shared<elephon::PhononStructure::Forces>();
	forces->initialize(displColl,
			res->get_electronic_structure_interface());

	elephon::PhononStructure::ForceConstantMatrix phi;
	phi.initialize( primitiveCell, supercell, displColl, primToSCConn, forces);

	BOOST_REQUIRE( phi.get_num_modes() == 3 );

	BOOST_REQUIRE( phi.get_num_R() == 2*2*2 );

	//Now check an explicit comparison for a comparable calculation with phonopy
	std::vector<double> referenceFC = {
			  3.64608298e+00,  -5.32878110e-17,  9.28043218e-17,
			 -1.49734409e-17,   3.64608298e+00, -6.41067535e-16,
			  1.70521532e-16,  -7.78792009e-16,  3.64608298e+00,

			 -8.65285861e-01,  -7.80371188e-18,  1.31237746e+00,
			 -2.15395352e-16,   1.93791099e-01,  1.51310358e-16,
			  1.31237746e+00,  -4.84520523e-17, -8.65285861e-01,

			  1.93791099e-01,  -2.06434132e-17, -2.15635582e-17,
			 -5.14092022e-18,  -8.65285861e-01, -1.31237746e+00,
			  1.32114617e-17,  -1.31237746e+00, -8.65285861e-01,

			 -8.65285861e-01,  -1.31237746e+00,  9.21050439e-18,
			 -1.31237746e+00,  -8.65285861e-01, -4.14649753e-17,
			  2.48136793e-16,   1.40682295e-16,  1.93791099e-01,

			 -8.65285861e-01,   1.31237746e+00, -5.36033273e-17,
			  1.31237746e+00,  -8.65285861e-01,  2.65530269e-17,
			 -3.06974720e-16,   2.34017976e-16,  1.93791099e-01,

			  1.93791099e-01,   3.45212010e-17, -4.23802400e-17,
			  5.65475126e-17,  -8.65285861e-01,  1.31237746e+00,
			 -5.30718169e-17,   1.31237746e+00, -8.65285861e-01,

			 -8.65285861e-01,   5.42598539e-17, -1.31237746e+00,
			  2.30589625e-16,   1.93791099e-01,  1.39660575e-16,
			 -1.31237746e+00,   1.91908168e-17, -8.65285861e-01,

			 -5.72518905e-01,   8.31101630e-18, -1.31398818e-17,
			  5.44854842e-18,  -5.72518905e-01,  1.05911250e-16,
			 -2.05835097e-17,   1.19537809e-16, -5.72518905e-01 };

	BOOST_CHECK_SMALL( phi(0, 0, 0, 0, 0) - referenceFC[  0] , 0.001 );
	BOOST_CHECK_SMALL( phi(0, 0, 0, 0, 1) - referenceFC[  1] , 0.001 );
	BOOST_CHECK_SMALL( phi(0, 0, 0, 0, 2) - referenceFC[  2] , 0.001 );
	BOOST_CHECK_SMALL( phi(0, 0, 0, 1, 0) - referenceFC[  3] , 0.001 );
	BOOST_CHECK_SMALL( phi(0, 0, 0, 1, 1) - referenceFC[  4] , 0.001 );
	BOOST_CHECK_SMALL( phi(0, 0, 0, 1, 2) - referenceFC[  5] , 0.001 );
	BOOST_CHECK_SMALL( phi(0, 0, 0, 2, 0) - referenceFC[  6] , 0.001 );
	BOOST_CHECK_SMALL( phi(0, 0, 0, 2, 1) - referenceFC[  7] , 0.001 );
	BOOST_CHECK_SMALL( phi(0, 0, 0, 2, 2) - referenceFC[  8] , 0.001 );

	BOOST_CHECK_SMALL( phi(1, 0, 0, 0, 0) - referenceFC[  9] , 0.001 );
	BOOST_CHECK_SMALL( phi(1, 0, 0, 0, 1) - referenceFC[ 10] , 0.001 );
	BOOST_CHECK_SMALL( phi(1, 0, 0, 0, 2) - referenceFC[ 11] , 0.001 );
	BOOST_CHECK_SMALL( phi(1, 0, 0, 1, 0) - referenceFC[ 12] , 0.001 );
	BOOST_CHECK_SMALL( phi(1, 0, 0, 1, 1) - referenceFC[ 13] , 0.001 );
	BOOST_CHECK_SMALL( phi(1, 0, 0, 1, 2) - referenceFC[ 14] , 0.001 );
	BOOST_CHECK_SMALL( phi(1, 0, 0, 2, 0) - referenceFC[ 15] , 0.001 );
	BOOST_CHECK_SMALL( phi(1, 0, 0, 2, 1) - referenceFC[ 16] , 0.001 );
	BOOST_CHECK_SMALL( phi(1, 0, 0, 2, 2) - referenceFC[ 17] , 0.001 );

	BOOST_CHECK_SMALL( phi(0, 1, 0, 0, 0) - referenceFC[ 18] , 0.001 );
	BOOST_CHECK_SMALL( phi(0, 1, 0, 0, 1) - referenceFC[ 19] , 0.001 );
	BOOST_CHECK_SMALL( phi(0, 1, 0, 0, 2) - referenceFC[ 20] , 0.001 );
	BOOST_CHECK_SMALL( phi(0, 1, 0, 1, 0) - referenceFC[ 21] , 0.001 );
	BOOST_CHECK_SMALL( phi(0, 1, 0, 1, 1) - referenceFC[ 22] , 0.001 );
	BOOST_CHECK_SMALL( phi(0, 1, 0, 1, 2) - referenceFC[ 23] , 0.001 );
	BOOST_CHECK_SMALL( phi(0, 1, 0, 2, 0) - referenceFC[ 24] , 0.001 );
	BOOST_CHECK_SMALL( phi(0, 1, 0, 2, 1) - referenceFC[ 25] , 0.001 );
	BOOST_CHECK_SMALL( phi(0, 1, 0, 2, 2) - referenceFC[ 26] , 0.001 );
}
BOOST_AUTO_TEST_SUITE_END()
