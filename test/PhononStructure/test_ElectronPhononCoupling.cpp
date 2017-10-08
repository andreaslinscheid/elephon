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
#include "fixtures/DataLoader.h"
#include "fixtures/scenarios.h"
#include <vector>

BOOST_AUTO_TEST_CASE( Gkkp_generate_regular_k_grid_q_zero )
{
	auto resHandl = elephon::test::fixtures::scenarios::load_Al_fcc_primitive_vasp_sc2x2x2();
	std::vector<double> kpts = { 0.0, 0.0, 0.0 };
	auto kppts = kpts;

	std::vector<int> bandsList = { 1, 2 };
	auto bandspList = bandsList;

	elephon::PhononStructure::ElectronPhononCoupling gkkp;
	gkkp.generate_gkkp_and_phonon(
			kpts, kppts,
			bandsList, bandspList,
			resHandl->get_phonon_obj(),
			resHandl->get_displacement_potential_obj(),
			resHandl->get_wfct_obj() );


}

