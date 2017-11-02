/*	This file test_BandStructureAnalysis.cpp is part of elephon.
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
 *  Created on: Sep 8, 2017
 *      Author: A. Linscheid
 */
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ElectronicStructure
#include <boost/test/unit_test.hpp>
#include "Auxillary/UnitConversion.h"
#include "IOMethods/VASPInterface.h"
#include "IOMethods/ResourceHandler.h"
#include "ElectronicStructure/BandStructureAnalysis.h"
#include "fixtures/MockStartup.h"
#include <vector>
#include <memory>
#include <fstream>
#include <limits>

std::vector<double>
read_mass_tens_file(boost::filesystem::path massTens)
{
	std::ifstream file(massTens.c_str(), std::ios::binary);
	if ( ! file.good() )
		throw std::runtime_error("cannot open mass tensor file");

	file.seekg(0, std::ios::beg);
	int size = file.tellg();
	file.seekg(0, std::ios::end);
	size = int(file.tellg()) - size;
	std::vector<char> buf(size);
	file.seekg(0, std::ios::beg);
	file.read(&buf[0], size);

	int numElem = size/sizeof(float);

	auto ptr = reinterpret_cast<float*>(buf.data());
	return std::vector<double>(ptr, ptr+numElem);
}

BOOST_AUTO_TEST_CASE( write_mass_tensor )
{
	using namespace elephon;
	test::fixtures::MockStartup ms;
	auto massTens = ms.get_data_for_testing_dir() / "phony" / "mass_tensor.dat";
	boost::filesystem::remove(massTens);

	LatticeStructure::RegularSymmetricGrid grid;
	LatticeStructure::Symmetry id;
	id.set_reciprocal_space_sym();
	grid.initialize(std::vector<int>{60, 60, 60}, 1e-6, std::vector<double>{0.0, 0.0, 0.0}, id);

	// create a cosine kx + cosine ky + cosine kz band structure
	auto g = grid.get_grid_dim();
	std::vector<double> bandData( g[0]*g[1]*g[2] );
	for (int k = 0 ; k < g[2]; ++k)
		for (int j = 0 ; j < g[1]; ++j)
			for (int i = 0 ; i < g[0]; ++i)
				bandData[i+g[0]*(j+g[1]*k)] = std::cos((2*M_PI/g[0])*i)+std::cos((2*M_PI/g[1])*j)+std::cos((2*M_PI/g[2])*k);

	ElectronicStructure::ElectronicBands bands;
	bands.initialize(1, 0.0, bandData, grid);

	ElectronicStructure::BandStructureAnalysis::write_mass_tensor_file(
			massTens.string(), bands, 0, std::vector<double>{-4.0, 4.0}, 0);

	BOOST_REQUIRE(boost::filesystem::exists(massTens));

	// this model has a local maximum at k = 0 0 0 and a local minima at k = -0.5 -0.5 -0.5
	std::ifstream file(massTens.c_str(), std::ios::binary);
	if ( ! file.good() )
		throw std::runtime_error("cannot open mass tensor file");

	file.seekg(0, std::ios::beg);
	int size = file.tellg();
	file.seekg(0, std::ios::end);
	size = int(file.tellg()) - size;
	std::vector<char> buf(size);
	file.seekg(0, std::ios::beg);
	file.read(&buf[0], size);

	// there is one maximum and one minimum in the first unit cell
	BOOST_REQUIRE_EQUAL( size/sizeof(float), 18*2 );

	int c = 0;
	std::vector<double> databuf(reinterpret_cast<float*>(buf.data()), reinterpret_cast<float*>(buf.data())+18*2);
	auto unity = Auxillary::units::INVERSE_EV_TIMES_A2_TO_ME;
	auto ptr = databuf.data();
	BOOST_CHECK_EQUAL(ptr[c++], 1); // maximum
	BOOST_CHECK_EQUAL(ptr[c++], 0); // band 1
	BOOST_CHECK_EQUAL(ptr[c++], 0.0); // kx = 0.0
	BOOST_CHECK_EQUAL(ptr[c++], 0.0); // ky = 0.0
	BOOST_CHECK_EQUAL(ptr[c++], 0.0); // kz = 0.0
	BOOST_CHECK_EQUAL(ptr[c++], 3.0); // energy maxium
	BOOST_CHECK_SMALL(ptr[c++]-(-unity), 1e-2); // eigenvalues, eigenvectors are degenerate and arbitrary
	c += 3;
	BOOST_CHECK_SMALL(ptr[c++]-(-unity), 1e-2); //
	c += 3;
	BOOST_CHECK_SMALL(ptr[c++]-(-unity), 1e-2); //
	c += 3;

	BOOST_CHECK_EQUAL(ptr[c++],-1); // minimum
	BOOST_CHECK_EQUAL(ptr[c++], 0); // band 1
	BOOST_CHECK_EQUAL(ptr[c++],-0.5); // kx =-0.5
	BOOST_CHECK_EQUAL(ptr[c++],-0.5); // ky =-0.5
	BOOST_CHECK_EQUAL(ptr[c++],-0.5); // kz =-0.5
	BOOST_CHECK_EQUAL(ptr[c++],-3.0); // energy maxium
	BOOST_CHECK_SMALL(ptr[c++]-unity, 1e-2); // eigenvalues, eigenvectors are degenerate and arbitrary
	c += 3;
	BOOST_CHECK_SMALL(ptr[c++]-unity, 1e-2); //
	c += 3;
	BOOST_CHECK_SMALL(ptr[c++]-unity, 1e-2); //
	c += 3;
	boost::filesystem::remove(massTens);
}

BOOST_AUTO_TEST_CASE( write_mass_tensor_skew_basis )
{
	using namespace elephon;
	test::fixtures::MockStartup ms;
	auto massTens = ms.get_data_for_testing_dir() / "phony" / "mass_tensor.dat";
	boost::filesystem::remove(massTens);

	std::vector<double>latticeMatrix{1.0/2.0, -1.0/2.0, 0, 1.0/2.0, 1.0/2.0, 0.0, 0.0, 0.0, 1.0};
	elephon::LatticeStructure::LatticeModule latticeMod( latticeMatrix );

	std::vector<int> g{64, 64, 64};
	int nG = g[0]*g[1]*g[2];
	int D = 3;
	std::vector<double> data(nG);

	// here we initialize data in the cubic cell. In terms of the lattice basis
	// this will of cause be tilted.
	std::vector<double> cubicVect(D);
	for ( int k = 0 ; k < g[2]; ++k)
		for ( int j = 0 ; j < g[1]; ++j)
			for ( int i = 0 ; i < g[0]; ++i)
			{
				cubicVect = {i/double(g[0]), j/double(g[1]), k/double(g[2])};
				latticeMod.reci_direct_to_cartesian_2pibya(cubicVect);
				int cnsq = i + g[0]*(j + g[1]*k);
				data[cnsq] = std::cos(cubicVect[0])+std::cos(cubicVect[1])
								+std::cos(cubicVect[2]);
			}

	// reinitialize such that the latticeMatrix above is the inverse matrix
	latticeMod.initialize(std::vector<double>{1/2.0, -(1/2.0), 0, 1/2.0, 1/2.0, 0.0, 0.0, 0.0, 1.0});

	LatticeStructure::RegularSymmetricGrid kgrid;
	LatticeStructure::Symmetry id;
	id.set_reciprocal_space_sym();
	kgrid.initialize(g, 1e-6, std::vector<double>{0.0, 0.0, 0.0}, id, latticeMod);

	ElectronicStructure::ElectronicBands bands;
	bands.initialize(1, 0.0, data, kgrid);

	ElectronicStructure::BandStructureAnalysis::write_mass_tensor_file(
			massTens.string(), bands, 0, std::vector<double>{-4.0, 4.0}, 0);

	BOOST_REQUIRE(boost::filesystem::exists(massTens));

	// this model has a local maximum at k = 0 0 0 and a local minima at k = -0.5 -0.5 -0.5
	std::ifstream file(massTens.c_str(), std::ios::binary);
	if ( ! file.good() )
		throw std::runtime_error("cannot open mass tensor file");

	file.seekg(0, std::ios::beg);
	int size = file.tellg();
	file.seekg(0, std::ios::end);
	size = int(file.tellg()) - size;
	std::vector<char> buf(size);
	file.seekg(0, std::ios::beg);
	file.read(&buf[0], size);

	// there are 2 maxima and 2 minima in the first unit cell
	BOOST_REQUIRE_EQUAL( size/sizeof(float), 18*4 );

	int c = 0;
	std::vector<double> databuf(reinterpret_cast<float*>(buf.data()), reinterpret_cast<float*>(buf.data())+18*4);
	auto ptr = databuf.data();
	auto unity = Auxillary::units::INVERSE_EV_TIMES_A2_TO_ME;
	BOOST_CHECK_EQUAL(ptr[c++], 1); // maximum
	BOOST_CHECK_EQUAL(ptr[c++], 0); // band 1
	BOOST_CHECK_EQUAL(ptr[c++], 0.0); // kx = 0.0
	BOOST_CHECK_EQUAL(ptr[c++], 0.0); // ky = 0.0
	BOOST_CHECK_EQUAL(ptr[c++], 0.0); // kz = 0.0
	BOOST_CHECK_EQUAL(ptr[c++], 3.0); // energy maxium
	BOOST_CHECK_SMALL(ptr[c++]-(-unity), 2e-2); // eigenvalues, eigenvectors are degenerate and arbitrary
	c += 3;
	BOOST_CHECK_SMALL(ptr[c++]-(-unity), 2e-2); //
	c += 3;
	BOOST_CHECK_SMALL(ptr[c++]-(-unity), 2e-2); //
	c += 3;

	BOOST_CHECK_EQUAL(ptr[c++], 1); // maximum
	BOOST_CHECK_EQUAL(ptr[c++], 0); // band 1
	BOOST_CHECK_EQUAL(ptr[c++],-0.5); // kx =-0.5
	BOOST_CHECK_EQUAL(ptr[c++],-0.5); // ky =-0.5
	BOOST_CHECK_EQUAL(ptr[c++], 0.0); // kz = 0.0
	BOOST_CHECK_EQUAL(ptr[c++], 3.0); // energy maxium
	BOOST_CHECK_SMALL(ptr[c++]-(-unity), 2e-2);
	c += 3;
	BOOST_CHECK_SMALL(ptr[c++]-(-unity), 2e-2);
	c += 3;
	BOOST_CHECK_SMALL(ptr[c++]-(-unity), 2e-2);
	c += 3;

	BOOST_CHECK_EQUAL(ptr[c++],-1); // minimum
	BOOST_CHECK_EQUAL(ptr[c++], 0); // band 1
	BOOST_CHECK_EQUAL(ptr[c++],-0.5); // kx =-0.5
	BOOST_CHECK_EQUAL(ptr[c++], 0.0); // ky = 0.0
	BOOST_CHECK_EQUAL(ptr[c++],-0.5); // kz =-0.5
	BOOST_CHECK_EQUAL(ptr[c++],-3.0); // energy maxium
	BOOST_CHECK_SMALL(ptr[c++]-unity, 2e-2);
	c += 3;
	BOOST_CHECK_SMALL(ptr[c++]-unity, 2e-2);
	c += 3;
	BOOST_CHECK_SMALL(ptr[c++]-unity, 2e-2);
	c += 3;

	BOOST_CHECK_EQUAL(ptr[c++],-1); // minimum
	BOOST_CHECK_EQUAL(ptr[c++], 0); // band 1
	BOOST_CHECK_EQUAL(ptr[c++], 0.0); // kx = 0.0
	BOOST_CHECK_EQUAL(ptr[c++],-0.5); // ky =-0.5
	BOOST_CHECK_EQUAL(ptr[c++],-0.5); // kz =-0.5
	BOOST_CHECK_EQUAL(ptr[c++],-3.0); // energy maxium
	BOOST_CHECK_SMALL(ptr[c++]-unity, 2e-2);
	c += 3;
	BOOST_CHECK_SMALL(ptr[c++]-unity, 2e-2);
	c += 3;
	BOOST_CHECK_SMALL(ptr[c++]-unity, 2e-2);
	c += 3;
	boost::filesystem::remove(massTens);
}

BOOST_AUTO_TEST_CASE( write_mass_tensor_fcc_Al_vasp )
{
	using namespace elephon;
	test::fixtures::MockStartup ms;
	auto testd = ms.get_data_for_testing_dir() / "Al" / "vasp" / "fcc_primitive";
	auto massTens = testd / "mass_tensor.dat";
	boost::filesystem::remove(massTens);
	std::string input = std::string()+
			"root_dir = "+testd.string()+"\n"
			"ewinbnd = -1.0 1.0\n"
			"f_mtens = "+massTens.string()+"\n";
	elephon::IOMethods::InputOptions opts;
	ms.simulate_elephon_input(
			(testd / "infile").string(),
			input,
			opts);

	auto loader = std::make_shared<elephon::IOMethods::VASPInterface>(opts);
	auto resources = std::make_shared<elephon::IOMethods::ResourceHandler>(loader);

	ElectronicStructure::BandStructureAnalysis::do_band_structure_analysis(resources);

	BOOST_REQUIRE(boost::filesystem::exists(massTens));

	boost::filesystem::remove(massTens);
}

BOOST_AUTO_TEST_CASE( Si_vasp_mass_tensor )
{
	using namespace elephon;
	test::fixtures::MockStartup ms;
	auto testd = ms.get_data_for_testing_dir() / "Si" / "vasp" / "primitive";
	auto massTens = testd / "mass_tensor.dat";
	boost::filesystem::remove(massTens);
	std::string input = std::string()+
			"root_dir = "+testd.string()+"\n"
			"ewinbnd = -4.0 4.0\n"
			"fftd = 64 64 64\n"
			"f_mtens = "+massTens.string()+"\n";
	elephon::IOMethods::InputOptions opts;
	ms.simulate_elephon_input(
			(testd / "infile").string(),
			input,
			opts);

	auto loader = std::make_shared<elephon::IOMethods::VASPInterface>(opts);
	auto resources = std::make_shared<elephon::IOMethods::ResourceHandler>(loader);

	ElectronicStructure::BandStructureAnalysis::do_band_structure_analysis(resources);

	BOOST_REQUIRE(boost::filesystem::exists(massTens));

	auto massTensorData = read_mass_tens_file(massTens);
	int numExtrema = massTensorData.size() / 18 ;
	std::vector< std::vector<double> > kPoints(numExtrema);
	std::vector< std::vector<double> > eigenvalues(numExtrema);

	// use a compare function that treats close doubles as equal
	auto cmp = [] (double a, double b) {
		if ( std::abs(a-b)>1e-3 )
			return a < b;
		return false;
		};
	std::multimap<double,int,decltype(cmp)> sortMinima(cmp);
	std::multimap<double,int,decltype(cmp)> sortMax(cmp);
	for (int ie = 0 ; ie < numExtrema ; ++ie )
	{
		kPoints[ie] = std::vector<double>(&massTensorData[18*ie+2], &massTensorData[18*ie+2]+3);
		eigenvalues[ie] = std::vector<double>{ massTensorData[18*ie+ 6],
											   massTensorData[18*ie+10],
											   massTensorData[18*ie+14] };
		double energy = massTensorData[18*ie+5];
		if ( massTensorData[18*ie+0] >= 0 )
		{
			sortMax.insert(std::make_pair(-energy, ie));
		}
		else
		{
			sortMinima.insert(std::make_pair(energy, ie));
		}
	}

	// find the valence band maxima and the conduction band minima
	std::vector<int> cbmin;
	double last = std::numeric_limits<double>::min();
	for ( auto min : sortMinima)
	{
		if ( min.first < 0 )
			continue;

		// min.first == last.first according to cmp
		if ( (!(cmp(last, min.first)) and (!cmp(min.first ,last))) or cbmin.empty() )
		{
			cbmin.push_back(min.second);
		}
		else
		{
			break;
		}
		last = min.first;
	}
	std::vector<int> vbmax;
	for ( auto max : sortMax)
	{
		if ( max.first < 0 ) // conduction bands; note negative energy convention
			continue;

		double energy = -max.first; // revert negative energy convention.

		// min.first == last.first according to cmp
		if ( (!(cmp(last, energy)) and (!cmp(energy, last))) or vbmax.empty() )
		{
			vbmax.push_back(max.second);
		}
		else
		{
			break;
		}
		last = energy;
	}

	for ( auto i : vbmax )
	{
		std::cout << "\nValence band max at k="<<
				kPoints[i][0] << ","<< kPoints[i][1] << ","<< kPoints[i][2] << " energy " << massTensorData[18*i+5]<<"\n";
		std::cout << "Mass tensor eigenvalues: "<<
				eigenvalues[i][0] << ","<< eigenvalues[i][1] << ","<< eigenvalues[i][2] <<"\n";
	}

	for ( auto i : cbmin )
	{
		std::cout << "\nConduction band min at k="<<
				kPoints[i][0] << ","<< kPoints[i][1] << ","<< kPoints[i][2] << " energy " << massTensorData[18*i+5]<<"\n";
		std::cout << "Mass tensor eigenvalues: "<<
				eigenvalues[i][0] << ","<< eigenvalues[i][1] << ","<< eigenvalues[i][2] <<"\n";
	}


	boost::filesystem::remove(massTens);
}
