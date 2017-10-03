/*	This file test_ElectronicBands.cpp is part of elephon.
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
 *  Created on: May 20, 2017
 *      Author: A. Linscheid
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Input_test
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "IOMethods/ReadVASPWaveFunction.h"
#include "IOMethods/ReadVASPSymmetries.h"
#include "IOMethods/ReadVASPPoscar.h"
#include "IOMethods/VASPInterface.h"
#include "Algorithms/FFTInterface.h"
#include <LatticeStructure/RegularSymmetricGrid.h>
#include "LatticeStructure/Symmetry.h"
#include "ElectronicStructure/ElectronicBands.h"
#include "fixtures/MockStartup.h"
#include "fixtures/DataLoader.h"
#include <vector>
#include <complex>
#include <cmath>

BOOST_AUTO_TEST_CASE( Bands_Symmetry_reconstruction )
{
	test::fixtures::MockStartup ms;
	auto testd = ms.get_data_for_testing_dir() / "FeSe" / "vasp" / "wfct";

	//Read data to construct both the band grid from both the symmetric and the non-symmetric input data.
	elephon::IOMethods::ReadVASPWaveFunction wfcSymRead;
	wfcSymRead.prepare_wavecar( (testd / "symmetric" / "WAVECAR").string() );

	elephon::IOMethods::ReadVASPWaveFunction wfcNoSymRead;
	wfcNoSymRead.prepare_wavecar( (testd / "no_symmetry" / "WAVECAR").string() );

	//First we make sure that the k point that irreducible k point have the same data in the reducible grid
	//The analogy is put for this hand crafted input data by hand!
	std::vector<int> irredNoSym = {0,1,2,6,7,12};
	assert( wfcSymRead.get_num_kpts() == 6 );
	BOOST_REQUIRE( wfcSymRead.get_num_bands() == wfcNoSymRead.get_num_bands());
	double diffIrreducibeWedge = 0;
	for ( int ik = 0 ; ik < wfcSymRead.get_num_kpts(); ++ik)
	{
		assert( std::abs(wfcSymRead.get_k_points()[ik*3+0]-wfcNoSymRead.get_k_points()[irredNoSym[ik]*3+0]) < 1e-6 );
		assert( std::abs(wfcSymRead.get_k_points()[ik*3+1]-wfcNoSymRead.get_k_points()[irredNoSym[ik]*3+1]) < 1e-6 );
		assert( std::abs(wfcSymRead.get_k_points()[ik*3+2]-wfcNoSymRead.get_k_points()[irredNoSym[ik]*3+2]) < 1e-6 );
		for ( int i = 0 ; i < wfcSymRead.get_num_bands(); ++i)
		{
			diffIrreducibeWedge += 1.0/wfcSymRead.get_num_bands()/3*std::abs(
					wfcSymRead.get_energies()[ik*wfcSymRead.get_num_bands()+i]
					-wfcNoSymRead.get_energies()[irredNoSym[ik]*wfcNoSymRead.get_num_bands()+i] );
		}
	}
	BOOST_REQUIRE( diffIrreducibeWedge < 1e-6 );

	std::string input = std::string()+
			"root_dir = "+(testd / "symmetric").string()+"\n";
	elephon::IOMethods::InputOptions opts;
	ms.simulate_elephon_input(
			(testd / "infile").string(),
			input,
			opts);

	auto loader = std::make_shared<elephon::IOMethods::VASPInterface>(opts);
	elephon::ElectronicStructure::ElectronicBands bands_sym;
	loader->read_band_structure((testd / "symmetric").string(), bands_sym);
	bands_sym.initialize( wfcSymRead.get_num_bands(), 0.0, wfcSymRead.get_energies(), bands_sym.get_grid() );

	elephon::LatticeStructure::Symmetry identity;
	identity.initialize(
			1e-6,
			std::vector<int>({1,0,0,0,1,0,0,0,1}),
			std::vector<double>({0,0,0}),
			bands_sym.get_grid().get_lattice(),
			false );
	identity.set_reciprocal_space_sym();

	elephon::LatticeStructure::RegularSymmetricGrid gridNoSym;
	gridNoSym.initialize(
			std::vector<int>({5,5,2}),
			1e-6,
			std::vector<double>({0.0,0.0,0.5}),
			identity,
			bands_sym.get_grid().get_lattice() );

	elephon::ElectronicStructure::ElectronicBands bands_nosym;
	bands_nosym.initialize( wfcNoSymRead.get_num_bands(), 0.0, wfcNoSymRead.get_energies(), gridNoSym );

	//Check that the reducible grids match
	BOOST_REQUIRE( bands_sym.get_grid().get_np_red() == bands_nosym.get_grid().get_np_red() );
	BOOST_REQUIRE( bands_sym.get_nBnd() == bands_nosym.get_nBnd() );

	std::vector<int> bndRequests( bands_sym.get_nBnd() );
	for ( int ibnd = 0 ; ibnd < bands_sym.get_nBnd();  ++ibnd )
		bndRequests[ibnd] = ibnd;
	std::vector<double> dataSym, dataNoSym;
	bands_sym.generate_reducible_data( bndRequests ,dataSym );
	bands_nosym.generate_reducible_data( bndRequests ,dataNoSym );

	double diff = 0;
	for ( int i = 0 ; i < bands_sym.get_grid().get_np_red();  ++i )
		for ( int j = 0 ; j < bands_sym.get_nBnd();  ++j )
			diff += 1.0/wfcSymRead.get_num_bands()/bands_sym.get_grid().get_np_red()*
					std::abs( dataSym[ i*bands_sym.get_nBnd()+j ]
							- dataNoSym[ i*bands_sym.get_nBnd()+j ] );
	BOOST_REQUIRE( diff < 1e-6 );
}

/** Test the generation of data on a reducible grid
 */
BOOST_AUTO_TEST_CASE( Bands_reducible_reconstruction)
{
	test::fixtures::DataLoader dl;
	auto bands = dl.create_symmetric_cosine_model({50, 50, 50}, {0.0, 0.0, 0.0});
	auto griddims = bands.get_grid().get_grid_dim();
	auto gridshift = bands.get_grid().get_grid_shift();
	int nB = bands.get_nBnd();
	assert(nB == 2);

	std::vector<double> bandDataReconstr;
	bands.generate_reducible_data(std::vector<int>{0, 1}, bandDataReconstr);
	BOOST_REQUIRE(bandDataReconstr.size() == nB*griddims[2]*griddims[1]*griddims[0]);
	double diff = 0.0;
	for ( int iz = 0 ; iz < griddims[2]; ++iz )
		for ( int iy = 0 ; iy < griddims[1]; ++iy )
			for ( int ix = 0 ; ix < griddims[0]; ++ix )
			{
				int cnsq = ix + griddims[0]*(iy + griddims[1]*iz);
				diff += std::abs(bandDataReconstr[cnsq*nB+0] - (std::cos((2*M_PI/griddims[0])*(ix+gridshift[0]))
								 +std::cos((2*M_PI/griddims[1])*(iy+gridshift[1]))));

				diff += std::abs(bandDataReconstr[cnsq*nB+1] - (std::cos((2*M_PI/griddims[2])*(iz+gridshift[2]))));
			}
	BOOST_CHECK_SMALL(diff/(nB*griddims[2]*griddims[1]*griddims[0]), 1e-6);
}

BOOST_AUTO_TEST_CASE( Bands_fft_interpolation )
{
	test::fixtures::DataLoader dl;
	auto bands = dl.create_symmetric_cosine_model({50, 50, 50}, {0.0, 0.0, 0.0});
	auto griddims = bands.get_grid().get_grid_dim();
	auto gridshift = bands.get_grid().get_grid_shift();
	int nB = bands.get_nBnd();
	assert(nB == 2);

	std::vector<int> griddimsNew{64, 64, 64};
	std::vector<double> gridshiftNew{0.5, 0.5, 0.5};
	bands.fft_interpolate(griddimsNew, gridshiftNew);
	auto const & kgrid = bands.get_grid();

	double diff = 0;
	for ( int iz = 0 ; iz < griddimsNew[2]; ++iz )
		for ( int iy = 0 ; iy < griddimsNew[1]; ++iy )
			for ( int ix = 0 ; ix < griddimsNew[0]; ++ix )
			{
				int cnsq = ix + griddimsNew[0]*(iy + griddimsNew[1]*iz);
				int ikir = kgrid.get_maps_red_to_irreducible()[cnsq];
				diff += std::abs(bands(ikir, 0) - std::cos((2*M_PI/griddimsNew[0])*(ix+gridshiftNew[0]))
								 -std::cos((2*M_PI/griddimsNew[1])*(iy+gridshiftNew[1])));
				diff += std::abs(bands(ikir, 1) -
									std::cos((2*M_PI/griddimsNew[2])*(iz+gridshiftNew[2])));
			}
	diff /= griddimsNew[0]*griddimsNew[1]*griddimsNew[2]*nB;

	BOOST_CHECK_SMALL(diff, 1e-6);
}

BOOST_AUTO_TEST_CASE( fft_interpol_MgB2_vasp )
{
	test::fixtures::MockStartup ms;
	auto testd = ms.get_data_for_testing_dir() / "MgB2" / "vasp" / "ldos";

	std::string input = std::string()+
			"root_dir = "+testd.string()+"\n";
	elephon::IOMethods::InputOptions opts;
	ms.simulate_elephon_input(
			(testd / "infile").string(),
			input,
			opts);

	auto loader = std::make_shared<elephon::IOMethods::VASPInterface>(opts);

	elephon::ElectronicStructure::ElectronicBands bands;
	loader->read_band_structure(testd.string(), bands);

	std::vector<int> gridd{32, 32, 32};
	std::vector<double> gshift{0.0, 0.0, 0.0}; // this is a hexagonal system - non-zero shift breaks the symmetry
	bands.fft_interpolate(gridd, gshift);

	BOOST_REQUIRE( bands.get_grid().get_grid_dim() == gridd );
}

BOOST_AUTO_TEST_CASE( fft_interpol_Al_vasp )
{
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

	elephon::ElectronicStructure::ElectronicBands bands;
	loader->read_band_structure(testd.string(), bands);

	std::vector<int> gridd{31, 31, 31};
	std::vector<double> gshift{0.0, 0.0, 0.0};
	bands.fft_interpolate(gridd, gshift);

	BOOST_REQUIRE( bands.get_grid().get_grid_dim() == gridd );
}
