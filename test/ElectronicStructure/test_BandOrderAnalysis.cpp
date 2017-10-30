/*	This file test_BandOrderAnalysis.cpp is part of elephon.
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
 *  Created on: Sep 12, 2017
 *      Author: A. Linscheid
 */
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ElectronicStructure
#include <boost/test/unit_test.hpp>
#include "ElectronicStructure/BandOrderAnalysis.h"
#include "fixtures/MockStartup.h"
#include "IOMethods/VASPInterface.h"
#include "IOMethods/ResourceHandler.h"
#include <vector>
#include <memory>
#include <iostream>

BOOST_AUTO_TEST_CASE( Al_fcc_vasp )
{
	exit(0);
	// the algorithm has to be redesigned.
	using namespace elephon;
	test::fixtures::MockStartup ms;
	auto testd = ms.get_data_for_testing_dir() / "Al" / "vasp" / "fcc_primitive";
	std::string input = std::string()+
			"root_dir = "+testd.string()+"\n";
	elephon::IOMethods::InputOptions opts;
	ms.simulate_elephon_input(
			(testd / "infile").string(),
			input,
			opts);

	auto loader = std::make_shared<elephon::IOMethods::VASPInterface>(opts);
	auto resource = std::make_shared<elephon::IOMethods::ResourceHandler>(loader);
	auto bands = resource->get_dense_electronic_bands_obj();
	auto wfcts = resource->get_wfct_obj();

	ElectronicStructure::BandOrderAnalysis analysis;
	analysis.compute_band_order_overlap(*bands, *wfcts);

	std::vector<double> regBndData;
	std::vector<int> bndIndices{0,1,2,3,4};
	bands->generate_reducible_data(bndIndices, regBndData);
	for ( int j = 0 ; j < 30; ++j)
	{
		for ( int i = 0 ; i < 30; ++i)
		{
			int ir = bands->get_grid().get_maps_red_to_irreducible()[i+30*j];
			std::cout << i << '\t' << j << '\t' <<
					regBndData[(i+30*j)*bndIndices.size()] << '\t' <<
				regBndData[(i+30*j)*bndIndices.size()+analysis(ir, 0)] << std::endl;
		}
	}
}
