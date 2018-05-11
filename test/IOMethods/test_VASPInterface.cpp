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
#include "Algorithms/helperfunctions.hpp"
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
		elephon::LatticeStructure::RegularBareGrid const & grid,
		elephon::Auxillary::alignedvector::DV const & regularGridPotential,
		elephon::AtomicSite::AtomSiteData const & radialPotential)
{
	std::ofstream file(filename.c_str());
	elephon::AtomicSite::RadialGrid const & rgrid = radialPotential.get_potential_data().get_radial_grid();

	std::vector<double> insideSphere;
	for (int i = 0 ; i < grid.get_grid_dim()[0]; ++i)
		for (int j = 0 ; j < grid.get_grid_dim()[1]; ++j)
		{
			std::vector<int> xyz({i,j,0});
			auto vector = grid.get_vector_direct(grid.get_xyz_to_reducible(xyz));
			grid.get_lattice().direct_to_cartesian_angstroem(vector);
			double r = std::sqrt(std::pow(vector[0]-rgrid.get_center()[0],2)+std::pow(vector[1]-rgrid.get_center()[1],2)
								+std::pow(vector[2]-rgrid.get_center()[2],2));
			if (r <= rgrid.get_range_of_definition())
			{
				insideSphere.insert(insideSphere.end(),{vector[0]+rgrid.get_center()[0], vector[1]+rgrid.get_center()[1],
						vector[2]+rgrid.get_center()[2]});
			}
		}
	elephon::Auxillary::alignedvector::ZV outputData;
	auto she = radialPotential.get_potential_data();
	radialPotential.get_frozen_core_data().add_core_hartree_potential(she.begin(), she.end());
	radialPotential.get_frozen_core_data().add_core_potential(she.begin(), she.end());
	she.interpolate(insideSphere, outputData);

	std::vector<double> ortho1, ortho2;
	elephon::Algorithms::helperfunctions::cross_prod(grid.get_lattice().get_lattice_vector(0),
			grid.get_lattice().get_lattice_vector(1),ortho1);
	elephon::Algorithms::helperfunctions::cross_prod(grid.get_lattice().get_lattice_vector(0),ortho1,ortho2);
	double ortho2Norm = std::sqrt(std::pow(ortho2[0],2)+std::pow(ortho2[1],2)+std::pow(ortho2[2],2));
	for (double & a : ortho2)
		a /= ortho2Norm;

	int c = 0;
	for (int i = 0 ; i < grid.get_grid_dim()[0]; ++i)
	{
		for (int j = 0 ; j < grid.get_grid_dim()[1]; ++j)
		{
			std::vector<int> xyz({i,j,0});
			const int cnsq = grid.get_xyz_to_reducible(xyz);
			std::vector<double> vec = grid.get_vector_direct(cnsq);
			grid.get_lattice().direct_to_cartesian_angstroem(vec);
			double r = std::sqrt(std::pow(vec[0]-rgrid.get_center()[0],2)+std::pow(vec[1]-rgrid.get_center()[1],2)
								+std::pow(vec[2]-rgrid.get_center()[2],2));
			if (r <= rgrid.get_range_of_definition()){
				file << vec[0]<< ' '<< vec[1]<< ' ' << vec[2]<< ' '<< std::real(outputData[c++]) << '\n';
			}else{
				file << vec[0]<< ' '<< vec[1]<< ' ' << vec[2]<< ' '<< std::real(regularGridPotential[cnsq]) << '\n';
			}
		}
	}
	file.close();
}

/**
 * Write a file to be plotted in 1D with gnuplot.
 *
 * @param[in] filename				Data file to be written
 * @param[in] dims					Dimension of the regularGridPotential
 * @param[in] regularGridPotential	Data on the regular grid
 * @param[in] radialPotential		Data on the radial grid
 */
void write_file_human_inspect_1D(std::string const & filename,
		elephon::LatticeStructure::RegularBareGrid const & grid,
		elephon::Auxillary::alignedvector::DV const & regularGridPotential,
		elephon::AtomicSite::AtomSiteData const & radialPotential)
{
	std::ofstream file(filename.c_str());
	std::vector<double> vec;
	std::vector<int> xyz;
	const int nSamples  = 60;
	std::vector<double> vectors(nSamples*3, 0.0);
	for (int ix = 0 ; ix < nSamples; ++ix)
		for (int i = 0 ; i < 3; ++i)
			vectors[ix*3+i] = double(ix+0.5)/nSamples;

	grid.get_lattice().direct_to_cartesian_angstroem(vectors);

	auto she = radialPotential.get_potential_data();
	radialPotential.get_frozen_core_data().add_core_hartree_potential(she.begin(), she.end());
	radialPotential.get_frozen_core_data().add_core_potential(she.begin(), she.end());
	elephon::Auxillary::alignedvector::ZV outputData;
	she.interpolate(vectors, outputData);
	for (int ix = 0 ; ix < nSamples; ++ix)
	{
		double Rad = radialPotential.get_frozen_core_data().get_radial_grid().get_range_of_definition();
		double magn = std::sqrt(std::pow(vectors[ix*3+0],2)+std::pow(vectors[ix*3+1],2)+std::pow(vectors[ix*3+2],2));
		if (magn <= Rad)
		{
			file << std::real(outputData[ix]) << ' ' << 0 << '\n';
		}
		else
		{
			int ixG = std::floor(double(ix)*grid.get_grid_dim()[0]/double(nSamples));
			int iyG = std::floor(double(ix)*grid.get_grid_dim()[1]/double(nSamples));
			int izG = std::floor(double(ix)*grid.get_grid_dim()[2]/double(nSamples));
			file << std::real(regularGridPotential[grid.get_xyz_to_reducible({ixG,iyG,izG})]) << ' ' << 1 << '\n';
		}
	}
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
	elephon::LatticeStructure::LatticeModule latmod;
	latmod.initialize({	-2.0197391656407362, 0.0000000000000000, 2.0197391656407362,
		     	 	 	 0.0000000000000000, 2.0197391656407362, 2.0197391656407362,
						-2.0197391656407362, 2.0197391656407362, 0.0000000000000000});
	elephon::LatticeStructure::RegularBareGrid grid(fftDims, false, 1e-6, {0.0, 0.0, 0.0}, latmod);

	double firstVal = -.12760769678E+03;
	double butLastVal = -.10585117804E+03;

	BOOST_CHECK_SMALL( std::abs(regularGridPotential[0] - firstVal), 1e-5 );
	double testVal = *(++regularGridPotential.rbegin());
	BOOST_CHECK_SMALL( std::abs( testVal - butLastVal), 1e-5 );

	BOOST_REQUIRE_EQUAL(radialPotential.size(),1u);

	double Al_frozen_core_charge = radialPotential[0].get_frozen_core_data().total_electronic_charge();
	BOOST_CHECK_CLOSE( Al_frozen_core_charge, 10.0, 1e-4 );

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

//	write_file_human_inspect((testd / "gnuplot_data.dat").string(),grid,regularGridPotential,radialPotential[0]);

//	write_file_human_inspect_1D((testd / "gnuplot_data.dat").string(),grid,regularGridPotential,radialPotential[0]);
}

BOOST_AUTO_TEST_SUITE_END()
