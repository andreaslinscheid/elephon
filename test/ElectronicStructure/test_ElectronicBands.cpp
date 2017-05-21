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
#include "LatticeStructure/RegularGrid.h"
#include "LatticeStructure/Symmetry.h"
#include "ElectronicStructure/ElectronicBands.h"
#include <vector>
#include <complex>
#include <cmath>

BOOST_AUTO_TEST_CASE( Bands_Symmetry_reconstruction )
{
	boost::filesystem::path dir = boost::filesystem::path(__FILE__).parent_path();

	//Read data to construct both the band grid from both the symmetric and the non-symmetric input data.
	elephon::IOMethods::ReadVASPWaveFunction wfcSymRead;
	wfcSymRead.prepare_wavecar( (dir / "FeSe_WAVECAR_sym.dat").string() );

	elephon::IOMethods::ReadVASPWaveFunction wfcNoSymRead;
	wfcNoSymRead.prepare_wavecar( (dir / "FeSe_WAVECAR_nosym.dat").string() );

	//First we make sure that the k point that irreducible k point have the same data in the reducible grid
	//The analogy is put for this hand crafted input data by hand!
	std::vector<int> irredNoSym = {0,1,5};
	assert( wfcSymRead.get_num_kpts() == 3 );
	BOOST_REQUIRE( wfcSymRead.get_num_bands() == wfcNoSymRead.get_num_bands());
	double diffIrreducibeWedge = 0;
	for ( int ik = 0 ; ik < 3; ++ik)
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

	elephon::IOMethods::ReadVASPSymmetries symread;
	symread.read_file( (dir / "FeSe_symmetries.dat").string() );

	elephon::LatticeStructure::Symmetry sym;
	sym.initialize(1e-6, symread.get_symmetries(), symread.get_fractionTranslations() );
	sym.set_reciprocal_space_sym();

	elephon::IOMethods::ReadVASPPoscar latRead;
	latRead.read_file( (dir / "FeSe_POSCAR.dat" ).string() );

	elephon::LatticeStructure::LatticeModule lattice;
	lattice.initialize( latRead.get_lattice_matrix() );

	elephon::LatticeStructure::RegularGrid grid;
	grid.initialize( 1e-6, std::vector<int>({4,4,1}), std::vector<double>({0.5,0.5,0.0}), sym, lattice );

	elephon::ElectronicStructure::ElectronicBands bands_sym;
	bands_sym.initialize( wfcSymRead.get_k_points(), wfcSymRead.get_num_bands(), wfcSymRead.get_energies(), grid );

	elephon::LatticeStructure::Symmetry identity;
	identity.initialize( 1e-6, std::vector<int>({1,0,0,0,1,0,0,0,1}), std::vector<double>({0,0,0}) );
	identity.set_reciprocal_space_sym();

	elephon::LatticeStructure::RegularGrid gridNoSym;
	gridNoSym.initialize( 1e-6, std::vector<int>({4,4,1}), std::vector<double>({0.5,0.5,0.0}), identity, lattice );

	elephon::ElectronicStructure::ElectronicBands bands_nosym;
	bands_nosym.initialize( wfcNoSymRead.get_k_points(), wfcNoSymRead.get_num_bands(), wfcNoSymRead.get_energies(), gridNoSym );

	//Check that the reducible grids match
	BOOST_REQUIRE( bands_sym.get_grid().get_np_red() == bands_nosym.get_grid().get_np_red() );
	BOOST_REQUIRE( bands_sym.get_nBnd() == bands_nosym.get_nBnd() );

	std::vector<int> bndRequests( bands_sym.get_nBnd() );
	for ( int ibnd = 0 ; ibnd < bands_sym.get_nBnd();  ++ibnd )
		bndRequests[ibnd] = ibnd;
	std::vector<double> dataSym, dataNoSym;
	bands_sym.generate_reducible_grid_bands( bndRequests ,dataSym );
	bands_nosym.generate_reducible_grid_bands( bndRequests ,dataNoSym );

	double diff = 0;
	for ( int i = 0 ; i < bands_sym.get_grid().get_np_red();  ++i )
		for ( int j = 0 ; j < bands_sym.get_nBnd();  ++j )
			diff += 1.0/wfcSymRead.get_num_bands()/bands_sym.get_grid().get_np_red()*
					std::abs( dataSym[ i*bands_sym.get_nBnd()+j ]
							- dataNoSym[ i*bands_sym.get_nBnd()+j ] );
	BOOST_REQUIRE( diff < 1e-6 );
}
