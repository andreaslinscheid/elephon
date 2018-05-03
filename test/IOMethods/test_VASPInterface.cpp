/*	This file test_VASPInterface.cpp is part of elephon.
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
 *  Created on: Apr 28, 2018
 *      Author: A. Linscheid
 */
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include <IOMethods/VASPInterface.h>
#include "fixtures/DataLoader.h"
#include "fixtures/MockStartup.h"
#include "Auxillary/AlignedVector.h"
#include <assert.h>
#include <vector>

BOOST_AUTO_TEST_SUITE( VASPInterface )

/**
 * Write a file to be plotted in 2D with gnuplot.
 *
 * @param[in] filename				Data file to be written
 * @param[in] dims					Dimension of the regularGridPotential
 * @param[in] regularGridPotential	Data on the regular grid
 * @param[in] radialPotential		Data on the radial grid
 */
void write_file_human_inspect(std::string const & filename,
		std::vector<int> const & dims,
		elephon::Auxillary::alignedvector::DV const & regularGridPotential,
		elephon::AtomicSite::AtomSiteData const & radialPotential)
{
	std::ofstream file(filename.c_str());
	elephon::AtomicSite::RadialGrid const & rgrid = radialPotential.get_potential_data().get_radial_grid();


	std::vector<double> insideSphere;
	for (int i = 0 ; i < dims[0]; ++i)
		for (int j = 0 ; j < dims[1]; ++j)
		{
			double x = double(i)/dims[0];
			double y = double(j)/dims[1];
			x -= std::floor(x+0.5);
			y -= std::floor(y+0.5);
			double r = std::sqrt(std::pow(x-rgrid.get_center()[0],2)+std::pow(y-rgrid.get_center()[1],2));
			if (r <= rgrid.get_range_of_definition()){
				insideSphere.insert(insideSphere.end(),{x+rgrid.get_center()[0], y+rgrid.get_center()[1], rgrid.get_center()[2]});
			}
		}
	elephon::Auxillary::alignedvector::ZV outputData;
	radialPotential.get_potential_data().interpolate(insideSphere, outputData);

	int c = 0;
	for (int i = 0 ; i < dims[0]; ++i)
	{
		for (int j = 0 ; j < dims[1]; ++j)
		{
			double x = double(i)/dims[0];
			double y = double(j)/dims[1];
			x -= std::floor(x+0.5);
			y -= std::floor(y+0.5);
			double r = std::sqrt(std::pow(x-rgrid.get_center()[0],2)+std::pow(y-rgrid.get_center()[1],2));
			if (r <= rgrid.get_range_of_definition()){
				file <<x<< ' '<< y<< ' '<< std::real(outputData[c++]) << '\n';
			}else{
				file <<x<< ' '<< y<< ' '<< std::real(regularGridPotential[i+dims[0]*j]) << '\n';
			}
		}
		file << '\n';
	}
	file.close();
}

BOOST_AUTO_TEST_CASE( Read_potcar )
{
	elephon::test::fixtures::MockStartup ms;
	auto testd = ms.get_data_for_testing_dir() / "Al" / "vasp" / "fcc_primitive" / "sc_4x4x4";
	elephon::test::fixtures::DataLoader dl;
	auto vaspLoader = dl.create_vasp_loader("");

	std::vector<int> fftDims;
	elephon::Auxillary::alignedvector::DV regularGridPotential;
	std::vector<elephon::AtomicSite::AtomSiteData> radialPotential;
	vaspLoader->read_electronic_potential(testd.string(), fftDims, regularGridPotential, radialPotential);

	BOOST_REQUIRE_EQUAL(fftDims.size(), 3);
	BOOST_CHECK_EQUAL(fftDims[0], 32);
	BOOST_CHECK_EQUAL(fftDims[1], 32);
	BOOST_CHECK_EQUAL(fftDims[2], 32);

	double firstVal = -.12760769678E+03;
	double butLastVal = -.10585117804E+03;

	BOOST_CHECK_SMALL( std::abs(regularGridPotential[0] - firstVal), 1e-5 );
	double testVal = *(++regularGridPotential.rbegin());
	BOOST_CHECK_SMALL( std::abs( testVal - butLastVal), 1e-5 );

	BOOST_REQUIRE_EQUAL(radialPotential.size(),1u);

	// it is difficult to check the radial potential for correctness. We employ one hint that
	// the transformations and reading was done correctly by confirming that even though we
	// internally use complex spherical harmonics in elephon, the potential is still real
	double imagPart = 0.0;
	const int numPts = 30;
	elephon::test::fixtures::DataLoader::RegularGridAtom atomGrid(numPts, radialPotential[0].get_potential_data().get_radial_grid());
	elephon::Auxillary::alignedvector::ZV outputData;
	radialPotential[0].get_potential_data().interpolate(atomGrid.get_atom_grid(), outputData);

	for (auto dataEle : outputData)
		imagPart += std::imag(dataEle);

	BOOST_CHECK_SMALL(imagPart, 1e-5);

	write_file_human_inspect((testd / "gnuplot_data.dat").string(),fftDims,regularGridPotential,radialPotential[0]);
}

BOOST_AUTO_TEST_SUITE_END()
