/*	This file test_DataRegularGrid.cpp is part of elephon.
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
 *  Created on: Nov 1, 2017
 *      Author: A. Linscheid
 */

#define BOOST_TEST_MODULE Input_test
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include <string>
#include "IOMethods/VASPInterface.h"
#include "IOMethods/ResourceHandler.h"
#include "IOMethods/Input.h"
#include "fixtures/MockStartup.h"

BOOST_AUTO_TEST_CASE( Si_vasp_electron_bands )
{
	elephon::test::fixtures::MockStartup ms;
	auto testd = ms.get_data_for_testing_dir() / "Si" / "vasp" / "primitive";
	auto fkpath = testd / "kpath.dat";
	auto fbands = testd / "bands.dat";
	boost::filesystem::remove(fkpath);
	boost::filesystem::remove(fbands);

	ms.write_kpath_file_fcc(fkpath.string());

	std::string input = std::string()+
			"root_dir = "+testd.string()+"\n"
			"fftd = 0 0 0\n"
			"ewinbnd = -5.0 5.0\n"
			"f_kpath = "+fkpath.string()+"\n"
			"f_bands = "+fbands.string()+"\n";
	elephon::IOMethods::InputOptions opts;
	ms.simulate_elephon_input(
			(testd / "infile").string(),
			input,
			opts);

	auto loader = std::make_shared<elephon::IOMethods::VASPInterface>(opts);
	auto resource = std::make_shared<elephon::IOMethods::ResourceHandler>(loader);

	auto bands = resource->get_dense_electronic_bands_obj();
	auto kpath = resource->get_k_path();
	elephon::Auxillary::alignedvector::DV bandsAlongPath;
	int numBandsInWindow;
	bands->interpolate_bands_along_path( kpath->get_k_points(),
			resource->get_optns().get_ewinbnd(),
			bandsAlongPath,
			numBandsInWindow,
			resource->get_interpol_reci_tetra_mesh_obj());

	auto gnuplotFile = fbands.string()+".gp";
	kpath->produce_gnuplot_script_stable_particle(
			gnuplotFile,
			fbands.string(),
			"E(k)",
			bandsAlongPath,
			numBandsInWindow,
			bands->interpret_range(resource->get_optns().get_ewinbnd()));

	boost::filesystem::remove(boost::filesystem::path(gnuplotFile));
	boost::filesystem::remove(fkpath);
	boost::filesystem::remove(fbands);
}

