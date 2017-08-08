/*	This file test_ElectronPhononCoupling.cpp is part of elephon.
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
 *  Created on: Jul 4, 2017
 *      Author: A. Linscheid
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE PhononStructure
#include <boost/test/unit_test.hpp>
#include "PhononStructure/ElectronPhononCoupling.h"
#include "ElectronicStructure/Wavefunctions.h"
#include "fixtures/MockStartup.h"
#include "fixtures/FixtureForceConstant.h"
#include "fixtures/DataLoader.h"
#include <vector>


elephon::ElectronicStructure::Wavefunctions
load_wfct_Al_fcc_primitive_vasp_sc2x2x2()
{
	test::fixtures::MockStartup ms;
	auto rootDir = ms.get_data_for_testing_dir() / "Al" / "vasp" / "fcc_primitive" ;
	auto phononDir = rootDir / "phonon";
	std::string content = std::string()+
			"scell=2 2 2\n"
			"root_dir="+rootDir.string()+"\n"
			"elphd="+phononDir.string()+"\n"
			"";
	test::fixtures::DataLoader dl;
	auto loader = dl.create_vasp_loader( content );

	elephon::ElectronicStructure::Wavefunctions wfcts;
	wfcts.initialize( 1e-6, rootDir.string(), loader);
	return wfcts;
}

BOOST_AUTO_TEST_CASE( Gkkp_generate_regular_k_grid_q_zero )
{
	test::fixtures::FixtureForceConstant ffc;
	auto fc = ffc.compute_fc_for_Al_gamma();

	std::vector<double> masses = {26.9815385};

	elephon::PhononStructure::Phonon ph;
	ph.initialize( fc, masses );

	std::vector<double> kpts = { 0.0, 0.0, 0.0 };
	auto kppts = kpts;

	std::vector<int> bandsList = { 1, 2 };
	auto bandspList = bandsList;

	test::fixtures::FixtureForceConstant ff;
	elephon::PhononStructure::DisplacementPotential dvscf
		= ff.build_displ_pot_Al_fcc_primitive_vasp_sc2x2x2();

	auto wfcts = load_wfct_Al_fcc_primitive_vasp_sc2x2x2();

	elephon::PhononStructure::ElectronPhononCoupling gkkp;
	gkkp.generate_gkkp( kpts, kppts, bandsList, bandspList, ph, dvscf, wfcts );


}

