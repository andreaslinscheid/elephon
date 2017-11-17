/*	This file test_LocalDerivative.cpp is part of elephon.
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
 *  Created on: Nov 2, 2017
 *      Author: A. Linscheid
 */

#define BOOST_TEST_MODULE localDerivative
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "IOMethods/VASPInterface.h"
#include "IOMethods/ResourceHandler.h"
#include "fixtures/MockStartup.h"
#include "Algorithms/LocalDerivatives.h"

class Functor
{

};

BOOST_AUTO_TEST_CASE( Analytic_simple )
{
	//todo move some of the tests for the dos and mass tensor here. They essentially test this method.
//	// load the silicone fcc lattice system
//	using namespace elephon;
//	test::fixtures::MockStartup ms;
//	auto testd = ms.get_data_for_testing_dir() / "Si" / "vasp" / "primitive";
//
//	std::string input = std::string()+
//			"root_dir = "+testd.string()+"\n";
//	elephon::IOMethods::InputOptions opts;
//	ms.simulate_elephon_input(
//			(testd / "infile").string(),
//			input,
//			opts);
//	auto loader = std::make_shared<elephon::IOMethods::VASPInterface>(opts);
//	auto resources = std::make_shared<elephon::IOMethods::ResourceHandler>(loader);
//
//	LatticeStructure::RegularBareGrid bgrid = resources->get_electronic_bands_obj()->get_grid().view_bare_grid();
//
//	compute_derivatives_sqr_polynom
}

