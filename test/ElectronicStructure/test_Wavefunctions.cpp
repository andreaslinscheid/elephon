/*	This file test_Wavefunctions.cpp is part of elephon.
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
 *  Created on: May 21, 2017
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
#include "LatticeStructure/RegularBareGrid.h"
#include "ElectronicStructure/Wavefunctions.h"
#include "fixtures/MockStartup.h"
#include "fixtures/DataLoader.h"
#include <vector>
#include <complex>
#include <cmath>
#include <set>

BOOST_AUTO_TEST_CASE( wavefunctions_partial_load )
{
	test::fixtures::MockStartup ms;
	auto rootDir = ms.get_data_for_testing_dir() / "Al" / "vasp" / "fcc_primitive" ;

	std::vector<double> kvec = {0.0, 0.0, 0.0};
	//Here we test the use case where we look up the cube around a point and compute those wavefunction
	//Wave functions are on a regular grid and in G (reciprocal) space
	//First, load the cubes of wave functions for linear interpolation

	std::string content = std::string("root_dir=")+rootDir.string()+"\n";
	test::fixtures::DataLoader dl;
	auto loader = dl.create_vasp_loader( content );

	elephon::ElectronicStructure::Wavefunctions wfcts;
	wfcts.initialize( 1e-6, rootDir.string(), loader);

	std::vector<int> kAToCube;
	std::vector<elephon::LatticeStructure::RegularBareGrid::GridCube> gridCubes;
	wfcts.get_k_grid().compute_grid_cubes_surrounding_nongrid_points(
			kvec,kAToCube,gridCubes);

	//Generate a list of all occurring reducible grid vectors.
	std::vector<int> npwPerKAllRedPts;
	std::vector<std::vector< std::complex<float> > > wfctAllRedPts;
	std::set<int> redIndicesSet;
	for ( auto c : gridCubes )
		redIndicesSet.insert( c.cornerIndices_.begin(),  c.cornerIndices_.end() );
	std::vector<int> redIndices(redIndicesSet.begin(), redIndicesSet.end());
	wfcts.generate_reducible_grid_wfcts(
			std::vector<int>{1,2}, //bands
			redIndices,
			wfctAllRedPts,
			npwPerKAllRedPts);
}

BOOST_AUTO_TEST_CASE( FeSe_Wfct_Symmetry_reconstruction )
{
	test::fixtures::MockStartup ms;
	auto testd = ms.get_data_for_testing_dir() / "FeSe" / "vasp" / "wfct";

	elephon::IOMethods::InputOptions opts;
	ms.simulate_elephon_input((testd/"infile").string(),"\n",opts);

	//introduce the VASP data loader
	std::shared_ptr<elephon::IOMethods::VASPInterface> loaderSym =
			std::make_shared< elephon::IOMethods::VASPInterface >(opts);

	std::shared_ptr<elephon::IOMethods::VASPInterface> loaderNoSym =
			std::make_shared< elephon::IOMethods::VASPInterface >(opts);

	elephon::ElectronicStructure::Wavefunctions wfct_sym;
	wfct_sym.initialize( 1e-6, (testd / "symmetric").string(), loaderSym );

	elephon::ElectronicStructure::Wavefunctions wfct_nosym;
	wfct_nosym.initialize( 1e-6, (testd / "no_symmetry").string(), loaderNoSym );

	elephon::LatticeStructure::LatticeModule lattice;
	elephon::LatticeStructure::Symmetry sym;
	std::vector<elephon::LatticeStructure::Atom> atoms;
	elephon::LatticeStructure::RegularSymmetricGrid kgrid;
	loaderSym->read_cell_paramters( (testd / "symmetric").string() ,1e-6,kgrid,lattice,atoms,sym);

	//compare band 0
	std::vector<int> bndInd = {0,5,11,15};
	std::vector<int> kptInd(kgrid.get_np_red());
	for (int ik = 0 ; ik < kgrid.get_np_red(); ++ik)
		kptInd[ik] = ik;

	std::vector<std::vector<std::complex<float>>> wfctnosym;
	std::vector<int> npPerKNoSym;
	wfct_nosym.generate_reducible_grid_wfcts(bndInd,kptInd,wfctnosym,npPerKNoSym);

	std::vector<std::vector<std::complex<float>>> wfctsym;
	std::vector<int> npPerKSym;
	wfct_sym.generate_reducible_grid_wfcts(bndInd,kptInd,wfctsym,npPerKSym);

	BOOST_REQUIRE( wfctnosym.size() == wfctsym.size() );
	BOOST_REQUIRE( npPerKSym == npPerKNoSym );

	int nB = int(bndInd.size());

	std::vector<int> fftMax = loaderSym->get_max_fft_dims();

	for (int ikir = 0 ; ikir < kgrid.get_np_irred(); ++ikir)
	{
		/*
		 *  We have selected a k point that has a full star so that each symmetry
		 * 	operation must generate an complete new wave function and at a given
		 * 	k point there is no degenerate subspace. This reduces ambiguity and allows
		 * 	for an easier comparison
		 */
		if ( ikir != 4)
			continue;

		auto symIrToRed = kgrid.get_maps_sym_irred_to_reducible()[ikir];
		auto irToRed = kgrid.get_maps_irreducible_to_reducible()[ikir];

		int npw = npPerKSym[ irToRed[0] ];

		auto kIrred = kgrid.get_vector_direct( irToRed[0] );
		std::vector<std::vector<int>> fftMapIrred;
		loaderNoSym->compute_fourier_map(kIrred, fftMapIrred, loaderNoSym->get_optns().get_gPrec());

		for (int istar = 0 ; istar < int(symIrToRed.size()) ; ++istar)
		{
			int ikRed = irToRed[istar];
			BOOST_REQUIRE( npw == npPerKNoSym[ikRed]);

			auto kRed = kgrid.get_vector_direct( ikRed );
			std::vector<std::vector<int>> fftMapRed;
			loaderNoSym->compute_fourier_map(kRed, fftMapRed, loaderNoSym->get_optns().get_gPrec());

			auto symOp = kgrid.get_symmetry().get_sym_op( symIrToRed[istar] );
			auto tau = kgrid.get_symmetry().get_fractional_translation( symIrToRed[istar] );

			/*
			 * First, we make sure that the loaded wavefunction in the non-symmetrized
			 * grid at a reducible k point is in fact the rotated wavefunction from the
			 * corresponding irreducible grid point.
			 */
			for ( int i = 0; i < nB; ++i )
			{
				auto irred_ptr = &wfctnosym[irToRed[0]][i*npw];
				std::vector< std::complex<float> > wfc1( irred_ptr, irred_ptr + npw);

				auto red_ptr = &wfctnosym[ikRed][i*npw];
				std::vector< std::complex<float> > wfc2( red_ptr, red_ptr + npw);

				std::map<int,int> nlm;
				for (int ipw = 0 ; ipw < npw ; ++ipw)
				{
					int csq =  fftMapRed[0][ipw*3] + fftMax[0]*(fftMapRed[0][ipw*3+1]+fftMax[1]*fftMapRed[0][ipw*3+2]);
					nlm.insert(std::make_pair(csq,ipw));
				}

				//choose the same phase convention
				int mpw = 0;
				float max = std::abs(wfc1[mpw]);
				for (int ipw = 0 ; ipw < npw ; ++ipw)
					if (std::abs(wfc1[ipw]) > max )
					{
						mpw = ipw;
						max = std::abs(wfc1[ipw]);
					}
				std::vector<double> Grot = {0,0,0};
				for ( int i = 0 ; i < 3 ; ++i)
					Grot[i] = symOp.ptgroup[i*3+0]*fftMapIrred[0][mpw*3]
							  + symOp.ptgroup[i*3+1]*fftMapIrred[0][mpw*3+1]
							 + symOp.ptgroup[i*3+2]*fftMapIrred[0][mpw*3+2];
				for ( int i = 0 ; i < 3 ; ++i)
					Grot[i] = Grot[i] < 0 ? Grot[i]+fftMax[i] : Grot[i];
				int csq = Grot[0] + fftMax[0]*(Grot[1]+fftMax[1]*Grot[2]);
				auto f = nlm.find(csq);

				int gx = fftMapIrred[0][mpw*3+0] < fftMax[0]/2 ? fftMapIrred[0][mpw*3+0]
																: fftMapIrred[0][mpw*3+0] - fftMax[0];
				int gy = fftMapIrred[0][mpw*3+1] < fftMax[1]/2 ? fftMapIrred[0][mpw*3+1]
																: fftMapIrred[0][mpw*3+1] - fftMax[1];
				int gz = fftMapIrred[0][mpw*3+2] < fftMax[2]/2 ? fftMapIrred[0][mpw*3+2]
																: fftMapIrred[0][mpw*3+2] - fftMax[2];
				std::complex<float> shiftPhase = std::exp(
						std::complex<float>(0,-2*M_PI*(tau[0]*gx+tau[1]*gy+tau[2]*gz)) );
				std::complex<float> phase = wfc1[mpw]/wfc2[f->second]/shiftPhase;
				for ( auto &c : wfc2 )
					c *= phase;

				//Make sure this is the same wave function
				float diff = 0;
				float norm1 = 0;
				float norm2 = 0;
				for (int ipw = 0 ; ipw < npw ; ++ipw)
				{
					auto Grot = kIrred;
					for ( int i = 0 ; i < 3 ; ++i)
						Grot[i] = symOp.ptgroup[i*3+0]*fftMapIrred[0][ipw*3]
								  + symOp.ptgroup[i*3+1]*fftMapIrred[0][ipw*3+1]
								 + symOp.ptgroup[i*3+2]*fftMapIrred[0][ipw*3+2];
					for ( int i = 0 ; i < 3 ; ++i)
						Grot[i] = Grot[i] < 0 ? Grot[i]+fftMax[i] : Grot[i];
					int csq = Grot[0] + fftMax[0]*(Grot[1]+fftMax[1]*Grot[2]);
					auto f = nlm.find(csq);
					assert ( f!= nlm.end() );

					int gx = fftMapIrred[0][ipw*3+0] < fftMax[0]/2 ? fftMapIrred[0][ipw*3+0]
																	: fftMapIrred[0][ipw*3+0] - fftMax[0];
					int gy = fftMapIrred[0][ipw*3+1] < fftMax[1]/2 ? fftMapIrred[0][ipw*3+1]
																	: fftMapIrred[0][ipw*3+1] - fftMax[1];
					int gz = fftMapIrred[0][ipw*3+2] < fftMax[2]/2 ? fftMapIrred[0][ipw*3+2]
																	: fftMapIrred[0][ipw*3+2] - fftMax[2];
					std::complex<float> shiftPhase =
							std::exp( std::complex<float>(0,-2*M_PI*(tau[0]*gx+tau[1]*gy+tau[2]*gz)) );
					diff += std::abs(wfc1[ipw]*shiftPhase-wfc2[f->second]);
					norm2 += std::abs(wfc2[f->second]);
					norm1 += std::abs(wfc1[ipw]);
				}
				BOOST_REQUIRE(  diff/(norm1+norm2) < 1e-3 );
			}

			/*
			 * Now we confirm that the data from the reducible (non symmetrized) calculation
			 *  and the calculation in the irreducible zone (symmetrized) that was rotated to
			 *  the reducible zone agrees.
			 */
			for ( int i = 0; i < nB; ++i )
			{
				auto ptrFromSym = &wfctsym[ikRed][i*npw];
				std::vector< std::complex<float> > wfc1( ptrFromSym, ptrFromSym + npw);

				auto ptrFromNoSym = &wfctnosym[ikRed][i*npw];
				std::vector< std::complex<float> > wfc2 ( ptrFromNoSym, ptrFromNoSym + npw);

				//Choose the same phase convention
				int mpw = 0;
				float max = std::abs(wfc1[mpw]);
				for (int ipw = 1 ; ipw < npw ; ++ipw)
					if (std::abs(wfc1[ipw]) > max )
					{
						mpw = ipw;
						max = std::abs(wfc1[ipw]);
					}

				//The difference is supposed to be just a phase, a.k.a |c|=1
				std::complex<float> phase = wfc1[mpw]/wfc2[mpw];
				BOOST_REQUIRE( std::abs(std::abs(phase)-1) < 1e-3 );
				for ( auto &c : wfc2 )
					c *= phase;

				//Make sure this is the same wave function
				float diff = 0;
				float norm1 = 0;
				float norm2 = 0;
				for (int ipw = 0 ; ipw < npw ; ++ipw)
				{
					diff += std::abs(wfc1[ipw]-wfc2[ipw]);
					norm2 += std::abs(wfc2[ipw]);
					norm1 += std::abs(wfc1[ipw]);
				}
				BOOST_REQUIRE(  diff/(norm1+norm2) < 1e-3 );
			}
		}
	}
}

BOOST_AUTO_TEST_CASE( Phony_VASP_Wfct_reconstruction )
{
	//We start by reading parameters of other files in this directory for lattice matrices and kpoints
	boost::filesystem::path phonyDir = boost::filesystem::path(__FILE__).parent_path()
											/ "../data_for_testing/phony/vasp_sym/";

	elephon::IOMethods::InputOptions noop;

	//introduce the VASP data loader
	std::shared_ptr<elephon::IOMethods::VASPInterface> loader =
			std::make_shared< elephon::IOMethods::VASPInterface >(noop);

	elephon::LatticeStructure::LatticeModule lattice;
	elephon::LatticeStructure::Symmetry sym;
	std::vector<elephon::LatticeStructure::Atom> atoms;
	elephon::LatticeStructure::RegularSymmetricGrid kgrid;
	loader->read_cell_paramters(phonyDir.string(),1e-6,kgrid,lattice,atoms,sym);
	sym.set_reciprocal_space_sym();

	//Now, we obtain some parameters from these function
	//Get the lattice basis and scale to their natural units
	auto B = lattice.get_reciprocal_latticeMatrix();
	for ( auto &bij : B )
		bij *= 2*M_PI/lattice.get_alat();
	const double eCut = 10;

	//Check the possible G vectors for plane wave coefficient

	//nGm should now be a lattice vector that is larger than the maximum plane
	//wave in any direction that is below the cutoff in energy.
	//Now we actually walk through them determine the maximum in each direction.
	int nGm = 100;
	const double energyConverionFactorVASP = 0.262465831;
	std::vector<int> fourierMax(3,0);
	std::vector<double> G(3);
	for ( int igx = -nGm; igx <= nGm; ++ igx)
		for ( int igy = -nGm; igy <= nGm; ++ igy)
			for ( int igz = -nGm; igz <= nGm; ++ igz)
			{
				for ( int i = 0 ; i < 3; ++i)
					G[i] = B[i*3+0]*igx+B[i*3+1]*igy+B[i*3+2]*igz;
				if ( (G[0]*G[0]+G[1]*G[1]+G[2]*G[2])/energyConverionFactorVASP < eCut )
				{
					fourierMax[0] = std::max(fourierMax[0],std::abs(igx));
					fourierMax[1] = std::max(fourierMax[1],std::abs(igy));
					fourierMax[2] = std::max(fourierMax[2],std::abs(igz));
				}
			}

	//We need to hold + and - direction (and zero). Also G+k can exceed this
	// number in each direction by up to one so we compute
	for ( auto & nGxi: fourierMax )
		nGxi = 2*(nGxi+1)+2;

	//We create a file that has the format of the VASP wavecar
	//but we control the input. VASP writes double precision (Fortran selected_real_kind(10))
	//or single precision (Fortran selected_real_kind(5)) dependent on the version of the code.
	typedef double VASPDP;
	typedef float VASPSP;

	std::ofstream file( (phonyDir / "WAVECAR" ).c_str() , std::ios_base::binary );
	if ( not file.good() )
		throw std::runtime_error( std::string("Error opening file ")+(phonyDir / "WAVECAR" ).string());

	//Parameters in the file
	const int nBnd = 1;
	const int nkpIrred = kgrid.get_np_irred();
	const int reclength = 100*sizeof(VASPDP);//In bytes

	//Compute the total number of records in the file as
	//the header plus per k point one record for bands ect and for each
	//band one record for the wavefunction data
	int numRecs = 2 + nkpIrred*(1+nBnd);
	std::vector<char> binaryBuffer( numRecs*reclength );
	VASPDP * buffAsFloat = reinterpret_cast<VASPDP * >( &binaryBuffer[0] );

	//Write header
	//First record
	buffAsFloat[0] = reclength;
	buffAsFloat[1] = 1; // #spin
	buffAsFloat[2] = 45200;//Version tag

	//second record
	buffAsFloat = reinterpret_cast<VASPDP * >( &binaryBuffer[reclength] );
	buffAsFloat[0] = nkpIrred;
	buffAsFloat[1] = nBnd;
	buffAsFloat[2] = eCut;
	for ( int i = 0 ; i < 3; ++i )
		for ( int j = 0 ; j < 3; ++j )
			buffAsFloat[3+i*3+j] = lattice.get_alat()*lattice.get_latticeMatrix()[i*3+j];

	//Header done. Now we create, write (and keep) the wavefunctions
	std::vector<std::vector<std::complex<float>>> wavefunctions(nkpIrred);
	std::vector<std::vector<int>> fftMaps(nkpIrred);
	std::vector<int> npwK(nkpIrred);

	for ( int ikir = 0 ; ikir < nkpIrred; ++ikir)
	{
		//Create the plane wave set at this k point
		int ikred = kgrid.get_maps_irreducible_to_reducible()[ikir][ kgrid.get_symmetry().get_identity_index() ];
		auto k = kgrid.get_vector_direct(ikred);
		auto kPlusG = k;
		fftMaps[ikir] = std::vector<int>(fourierMax[0]*fourierMax[1]*fourierMax[2]*3);
		int ng = 0;
		for ( int igz = 0 ; igz < fourierMax[2]; ++igz)
		{
			int igzf = igz < fourierMax[2]/2 ? igz : igz - fourierMax[2];
			for ( int igy = 0 ; igy < fourierMax[1]; ++igy)
			{
				int igyf = igy < fourierMax[1]/2 ? igy : igy - fourierMax[1];
				for ( int igx = 0 ; igx < fourierMax[0]; ++igx)
				{
					int igxf = igx < fourierMax[0]/2 ? igx : igx - fourierMax[0];

					for ( int xi = 0 ; xi < 3; ++xi)
						kPlusG[xi] = B[xi*3+0]*(k[0]+igxf) +B[xi*3+1]*(k[1]+igyf) +B[xi*3+2]*(k[2]+igzf);

					if ( (kPlusG[0]*kPlusG[0] + kPlusG[1]*kPlusG[1]
						  + kPlusG[2]*kPlusG[2])/energyConverionFactorVASP < eCut )
					{
						fftMaps[ikir][ng*3+0] = igx;
						fftMaps[ikir][ng*3+1] = igy;
						fftMaps[ikir][ng*3+2] = igz;
						ng++;
					}
				}
			}
		}
		npwK[ikir] = ng;
		fftMaps[ikir].resize(ng);

		wavefunctions[ikir].resize(ng*nBnd);
		for (int ibnd = 0 ; ibnd < nBnd; ++ibnd)
			for ( int ipw = 0 ; ipw < ng; ++ipw)
				wavefunctions[ikir][ibnd*ng+ipw] = std::complex<float>(ipw,ibnd);

		buffAsFloat = reinterpret_cast<VASPDP * >( &binaryBuffer[(2+ikir*(1+nBnd))*reclength] );
		buffAsFloat[0] = ng;
		buffAsFloat[1] = k[0];
		buffAsFloat[2] = k[1];
		buffAsFloat[3] = k[2];
		for (int ibnd = 0 ; ibnd < nBnd; ++ibnd)
			buffAsFloat[4+ibnd] = ibnd; // We set the energies to the band index

		for (int ibnd = 0 ; ibnd < nBnd; ++ibnd)
		{
			VASPSP * buffAsSingle = reinterpret_cast<VASPSP * >( &binaryBuffer[(2+ikir*(1+nBnd)+1+ibnd)*reclength] );
			for ( int ipw = 0 ; ipw < ng; ++ipw)
			{
				buffAsSingle[ipw*2+0] = std::real(wavefunctions[ikir][ibnd*ng+ipw]);
				buffAsSingle[ipw*2+1] = std::imag(wavefunctions[ikir][ibnd*ng+ipw]);
			}
		}
	}
	file.write(binaryBuffer.data(),binaryBuffer.size());
	file.close();

	//Read in the file once more and check if the write/read works
	std::vector<int> kpts( nkpIrred );
	for ( int ik = 0 ; ik < nkpIrred; ++ik)
		kpts[ik] = ik;
	std::vector<int> bandIndices( nBnd );
	for ( int ibnd = 0 ; ibnd < nBnd; ++ibnd)
		bandIndices[ibnd] = ibnd;
	std::vector<std::vector< std::complex<float> >> wfcts;
	std::vector<int> npwPerK;
	loader->read_wavefunctions( phonyDir.string(), kpts, bandIndices, wfcts, npwPerK );

	//Compare if we have read in the same thing we wrote to disk
	BOOST_REQUIRE( npwPerK == npwK);

	for ( int ik = 0 ; ik < nkpIrred; ++ik)
	{
		float diff = 0;
		for (int ibnd = 0 ; ibnd < nBnd; ++ibnd)
			for (int ipw = 0 ; ipw < nBnd; ++ipw)
				diff += std::abs(wavefunctions[ik][ibnd*npwK[ik]+ipw]-wfcts[ik][ibnd*npwK[ik]+ipw]);
		BOOST_REQUIRE( diff < 1e-5);
	}

	//Check if we have the irreducible wave function data as in our example
	auto compare_float = [] (float a, float b){
		return std::fabs(a-b) < 1e-5;
	};
	BOOST_REQUIRE( compare_float(std::real(wfcts[0][0]),0.0) && compare_float(std::imag(wfcts[0][0]),0.0) );
	BOOST_REQUIRE( compare_float(std::real(wfcts[0][1]),1.0) && compare_float(std::imag(wfcts[0][0]),0.0) );
	BOOST_REQUIRE( compare_float(std::real(wfcts[0][2]),2.0) && compare_float(std::imag(wfcts[0][0]),0.0) );
	BOOST_REQUIRE( compare_float(std::real(wfcts[0][3]),3.0) && compare_float(std::imag(wfcts[0][0]),0.0) );
	BOOST_REQUIRE( compare_float(std::real(wfcts[0][4]),4.0) && compare_float(std::imag(wfcts[0][0]),0.0) );

	//Seems so ... lets check the rotation
	elephon::ElectronicStructure::Wavefunctions phonyWave;
	phonyWave.initialize(1e-6,phonyDir.string(),loader);

	std::vector<int> reducibleIndices(kgrid.get_np_red());
	for (int i = 0 ; i < kgrid.get_np_red(); ++i)
		reducibleIndices[i] = i;
	std::vector<std::vector<std::complex<float>>> wfctsPerK;
	phonyWave.generate_reducible_grid_wfcts(bandIndices,reducibleIndices,wfctsPerK,npwPerK);

	for ( int ikir = 0 ; ikir < nkpIrred; ++ikir)
	{
		auto star = kgrid.get_maps_sym_irred_to_reducible();
		auto irredToRed = kgrid.get_maps_irreducible_to_reducible();
		for ( int is = 0 ; is < int(star.size()); ++is)
		{
			int isym = star[ikir][is];

			//Stars must have the same number of G vectors
			int ikreducible = kgrid.get_maps_irreducible_to_reducible()[ikir][is];
			BOOST_REQUIRE( npwPerK[ikreducible] == npwK[ikir] );

			//Check our particular example
			if ( (ikir == 0) and (isym == 3) )
			{
				BOOST_REQUIRE( compare_float(std::real(wfctsPerK[ikreducible][0]),0.0)
						&& compare_float(std::imag(wfctsPerK[ikreducible][0]),0.0) );
				BOOST_REQUIRE( compare_float(std::real(wfctsPerK[ikreducible][1]),2.0)
						&& compare_float(std::imag(wfctsPerK[ikreducible][1]),0.0) );
				BOOST_REQUIRE( compare_float(std::real(wfctsPerK[ikreducible][2]),1.0)
						&& compare_float(std::imag(wfctsPerK[ikreducible][2]),0.0) );
				BOOST_REQUIRE( compare_float(std::real(wfctsPerK[ikreducible][3]),4.0)
						&& compare_float(std::imag(wfctsPerK[ikreducible][3]),0.0) );
				BOOST_REQUIRE( compare_float(std::real(wfctsPerK[ikreducible][4]),3.0)
						&& compare_float(std::imag(wfctsPerK[ikreducible][4]),0.0) );
			}
		}
	}
}

BOOST_AUTO_TEST_CASE( MgB2_vasp_wfct_arbitray_kpts )
{
	test::fixtures::MockStartup ms;
	auto testd = ms.get_data_for_testing_dir() / "MgB2" / "vasp" / "ldos";
	auto outfile = testd / "ldos.dat";

	std::string input = std::string()+
			"root_dir = "+testd.string()+"\n";
	elephon::IOMethods::InputOptions opts;
	ms.simulate_elephon_input(
			(testd / "infile").string(),
			input,
			opts);

	auto loader = std::make_shared<elephon::IOMethods::VASPInterface>(opts);

	elephon::ElectronicStructure::Wavefunctions wfcts;
	wfcts.initialize(
			loader->get_optns().get_gPrec(),
			testd.string(),
			loader);

	std::vector<std::vector<std::complex<float>>> wfctsArbK;
	std::vector<std::vector<int>> fftMap;
	std::vector<double> k{0.0, 0.0, 0.0};
	std::vector<int> bands{0,1};
	wfcts.generate_wfcts_at_arbitray_kp(
			k,
			bands,
			wfctsArbK,
			fftMap);


}
