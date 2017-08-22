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
#include "IOMethods/WriteVASPWaveFunctions.h"
#include "IOMethods/ReadVASPSymmetries.h"
#include "IOMethods/ReadVASPPoscar.h"
#include "IOMethods/VASPInterface.h"
#include "LatticeStructure/RegularBareGrid.h"
#include "ElectronicStructure/Wavefunctions.h"
#include "fixtures/MockStartup.h"
#include "fixtures/DataLoader.h"
#include "Algorithms/FFTInterface.h"
#include "IOMethods/WriteVASPRealSpaceData.h"
#include <vector>
#include <cmath>
#include <complex>
#include <set>
#include <stdlib.h>
#include <time.h>

elephon::ElectronicStructure::Wavefunctions get_MgB2_vasp_wfct()
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
			testd.string(),
			loader);
	return wfcts;
}

void write_MgB2_vasp_chg(
		std::vector<int> const & chargeDim,
		std::vector<double> const & chgData)
{
	test::fixtures::DataLoader dl;
	test::fixtures::MockStartup ms;
	auto testd = ms.get_data_for_testing_dir() / "MgB2" / "vasp" / "ldos";
	std::string input = std::string()+
			"root_dir = "+testd.string()+"\n";
	auto uc = dl.load_unit_cell(input);
	elephon::IOMethods::WriteVASPRealSpaceData writer;
	writer.write_file(
			(testd / "chgGam.dat").string(),
			"charge due to wfcts at gamma",
			chargeDim,
			uc,
			chgData,
			false,
			true);
}

void write_Al_vasp_chg(
		std::vector<int> const & chargeDim,
		std::vector<double> const & chgData)
{
	test::fixtures::DataLoader dl;
	test::fixtures::MockStartup ms;
	auto testd = ms.get_data_for_testing_dir() / "Al" / "vasp" / "fcc_primitive" ;
	std::string input = std::string()+
			"root_dir = "+testd.string()+"\n";
	auto uc = dl.load_unit_cell(input);
	elephon::IOMethods::WriteVASPRealSpaceData writer;
	writer.write_file(
			(testd / "output"/ "LOCPOT").string(),
			"charge due to wfcts at gamma",
			chargeDim,
			uc,
			chgData,
			false,
			true);
}


double comp_symm_distortion(
		std::vector<double> const & tobechecked, //must be x-major
		std::vector<int> const & grid,
		elephon::LatticeStructure::Symmetry const & symmetry)
{
	assert( not symmetry.is_reci() );
	double result = 0.0;
	double integral = 0.0;
	for ( int k = 0; k < grid[2]; ++k)
		for ( int j = 0; j < grid[1]; ++j)
			for ( int i = 0; i < grid[0]; ++i)
			{
				double val = tobechecked[i+grid[0]*(j+grid[1]*k)];
				std::vector<double> v = {double(i)/grid[0], double(j)/grid[1], double(k)/grid[2]};
				for ( int is = 0; is < symmetry.get_num_symmetries(); ++is)
				{
					symmetry.apply(is, v, true);
					int ii = std::floor((v[0]-std::floor(v[0]))*grid[0]+0.5);
					int jj = std::floor((v[1]-std::floor(v[1]))*grid[1]+0.5);
					int kk = std::floor((v[2]-std::floor(v[2]))*grid[2]+0.5);
					int cnsq = ii+grid[0]*(jj+grid[1]*kk);
					assert((cnsq >=0) && (cnsq < tobechecked.size()));
					result += std::abs(val-tobechecked[cnsq]);
					integral += tobechecked[cnsq];
				}
			}
	return result/integral;
}

elephon::ElectronicStructure::Wavefunctions
load_Al_fcc_vasp_wfcts()
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
	wfcts.initialize(rootDir.string(), loader);
	return wfcts;
}

BOOST_AUTO_TEST_CASE( Al_vasp_wfct_interpol_star )
{
	auto wfcts = load_Al_fcc_vasp_wfcts();

	std::vector<int> chargeDim = {32, 32, 32};
	auto kg = wfcts.get_k_grid();
	auto rsSym = kg.get_symmetry();
	rsSym.set_reciprocal_space_sym(false);

	//Generate a random k vector and its star.
//	srand (time(NULL));
//	std::vector<double> kb{	float(rand()) / float(RAND_MAX),
//							float(rand()) / float(RAND_MAX),
//							float(rand()) / float(RAND_MAX) };
	// Since we do not have an upper bound for the allowed symmetry deviation testing a random k vector
	// is not sensible ... unfortunately.
	std::vector<double> kb{0.0967329, -0.212773, 0.274938};
	for ( auto & kxi : kb )
		kxi -= std::floor(kxi+0.5);
	auto genStar = kg.get_symmetry().star_operations( kb );
	std::vector<double> k(3+genStar.size()*3);
	std::copy(&kb[0], &kb[0]+3, &k[0]);
	for ( int istar = 0; istar < genStar.size(); ++istar )
	{
		auto kr = kb;
		genStar[istar].apply(kr);
		std::copy(&kr[0], &kr[0]+3, &k[3+istar*3]);
	}

	std::vector<int> bands{0};
	std::vector<std::vector<std::complex<float>>> wfctsArbK;
	std::vector<std::vector<int>> fftMap;
	wfcts.generate_wfcts_at_arbitray_kp(
			k,
			bands,
			wfctsArbK,
			fftMap);

	std::vector<std::complex<float>> wfctsGammaFullGrid;
	elephon::Algorithms::FFTInterface fft;

	std::vector<double> chggam(chargeDim[0]*chargeDim[1]*chargeDim[2], 0.0);
	for ( int ik = 0; ik < k.size()/3; ++ik)
	{
		fft.fft_sparse_data(
				fftMap[ik],
				wfcts.get_max_fft_dims(),
				wfctsArbK[ik],
				bands.size(),
				-1,
				wfctsGammaFullGrid,
				chargeDim,
				false,
				k.size()/3);

		for ( int i = 0; i < chggam.size(); ++i )
			chggam[i] += std::pow(std::real(wfctsGammaFullGrid[i]),2)+std::pow(std::imag(wfctsGammaFullGrid[i]),2);
	}
	write_Al_vasp_chg(chargeDim, chggam);
	auto r = comp_symm_distortion(chggam, chargeDim, rsSym);
	// Strictly speaking this method can slightly break the symmetry.
	// This is because the corner points of a cube are the symmetric wavefunctions. If, e.g., a point
	// is rotated to the same cube, the two wavefunctions are not correctly symmetry related.
	// Since, it is not obvious (to me) if or how one can establish an upper bound,
	// I am putting 5%, but actually this test is questionable since it may fail even
	// though the code is correct. NOTE: fixed k vector to keep this test.
	std::cout << "Testing k vector = " << kb[0] << ' ' << kb[1] << ' '<< kb[2]
			  << " and its star; Symmetry deviation = " << r*100 << "%"<< std::endl;
	BOOST_CHECK_SMALL(r, 5e-2);
}

BOOST_AUTO_TEST_CASE( wavefunctions_partial_load )
{
	//Here we test the use case where we look up the cube around a point and compute those wavefunction
	//Wave functions are on a regular grid and in G (reciprocal) space
	//First, load the cubes of wave functions for linear interpolation
	auto wfcts = load_Al_fcc_vasp_wfcts();

	std::vector<double> kvec = {0.0, 0.0, 0.0};
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
	wfct_sym.initialize((testd / "symmetric").string(), loaderSym );

	elephon::ElectronicStructure::Wavefunctions wfct_nosym;
	wfct_nosym.initialize((testd / "no_symmetry").string(), loaderNoSym );

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
	test::fixtures::MockStartup ms;
	boost::filesystem::path phonyDir = ms.get_data_for_testing_dir() / "phony" / "vasp_sym";

	elephon::IOMethods::InputOptions opts;
	ms.simulate_elephon_input((phonyDir/"infile").string(),"\n",opts);

	//introduce the VASP data loader
	std::shared_ptr<elephon::IOMethods::VASPInterface> loader =
			std::make_shared< elephon::IOMethods::VASPInterface >(opts);

	elephon::LatticeStructure::LatticeModule lattice;
	elephon::LatticeStructure::Symmetry sym;
	std::vector<elephon::LatticeStructure::Atom> atoms;
	elephon::LatticeStructure::RegularSymmetricGrid kgrid;
	loader->read_cell_paramters(phonyDir.string(),1e-6,kgrid,lattice,atoms,sym);

	const int nBnd = 1;
	const double eCut = 10;
	std::vector<double> bandData(kgrid.get_np_irred(), 0.0);
	elephon::ElectronicStructure::ElectronicBands bands;
	bands.initialize(nBnd, 0.0, bandData, kgrid);

	elephon::IOMethods::ReadVASPWaveFunction reader;
	std::vector<int> mFFT;
	reader.compute_fourier_max(eCut, lattice, mFFT);
	std::vector<std::vector<std::complex<float>>> wavefunctions(kgrid.get_np_irred());
	std::vector<std::vector<int>> fftMap;
	std::vector<double> kvectors(3*kgrid.get_np_irred());
	for (int ikir = 0 ; ikir < kgrid.get_np_irred(); ++ikir)
	{
		auto k = kgrid.get_vector_direct(kgrid.get_maps_irreducible_to_reducible()[ikir][kgrid.get_symmetry().get_identity_index()]);
		kvectors[ikir*3+0] = k[0];
		kvectors[ikir*3+1] = k[1];
		kvectors[ikir*3+2] = k[2];
	}
	reader.compute_fourier_map(
			kvectors,
			fftMap,
			kgrid.get_grid_prec(),
			1,
			mFFT,
			eCut,
			lattice);

	std::vector<int> npwK(kgrid.get_np_irred());
	for (int ikir = 0 ; ikir < kgrid.get_np_irred(); ++ikir)
	{
		npwK[ikir] = fftMap[ikir].size()/3;
		wavefunctions[ikir].resize(npwK[ikir]*nBnd);
		for (int ibnd = 0 ; ibnd < nBnd; ++ibnd)
			for ( int ipw = 0 ; ipw < npwK[ikir]; ++ipw)
				wavefunctions[ikir][ibnd*npwK[ikir] + ipw] = std::complex<float>(ipw,ibnd);
	}

	elephon::IOMethods::write_VASP_wavefunctions(
			(phonyDir / "WAVECAR" ).string(),
			eCut,
			lattice,
			bands,
			wavefunctions,
			fftMap);

	//Read in the file once more and check if the write/read works
	std::vector<int> kpts( kgrid.get_np_irred() );
	for ( int ik = 0 ; ik < kgrid.get_np_irred(); ++ik)
		kpts[ik] = ik;
	std::vector<int> bandIndices( nBnd );
	for ( int ibnd = 0 ; ibnd < nBnd; ++ibnd)
		bandIndices[ibnd] = ibnd;
	std::vector<std::vector< std::complex<float> >> wfcts;
	std::vector<int> npwPerK;
	loader->read_wavefunctions( phonyDir.string(), kpts, bandIndices, wfcts, npwPerK );

	//Compare if we have read in the same thing we wrote to disk
	BOOST_REQUIRE( npwPerK == npwK);

	for ( int ik = 0 ; ik < kgrid.get_np_irred(); ++ik)
	{
		float diff = 0;
		for (int ibnd = 0 ; ibnd < nBnd; ++ibnd)
			for (int ipw = 0 ; ipw < nBnd; ++ipw)
				diff += std::abs(wavefunctions[ik][ibnd*npwK[ik]+ipw]-wfcts[ik][ibnd*npwK[ik]+ipw]);
		BOOST_CHECK_SMALL( diff, float(1e-5));
	}

	//Check if we have the irreducible wave function data as in our example
	auto compare_float = [] (float a, float b){
		return std::fabs(a-b) < 1e-5;
	};
	BOOST_CHECK( compare_float(std::real(wfcts[0][0]),0.0) && compare_float(std::imag(wfcts[0][0]),0.0) );
	BOOST_CHECK( compare_float(std::real(wfcts[0][1]),1.0) && compare_float(std::imag(wfcts[0][0]),0.0) );
	BOOST_CHECK( compare_float(std::real(wfcts[0][2]),2.0) && compare_float(std::imag(wfcts[0][0]),0.0) );
	BOOST_CHECK( compare_float(std::real(wfcts[0][3]),3.0) && compare_float(std::imag(wfcts[0][0]),0.0) );
	BOOST_CHECK( compare_float(std::real(wfcts[0][4]),4.0) && compare_float(std::imag(wfcts[0][0]),0.0) );

	//Seems so ... lets check the rotation
	elephon::ElectronicStructure::Wavefunctions phonyWave;
	phonyWave.initialize(phonyDir.string(), loader);

	std::vector<int> reducibleIndices(kgrid.get_np_red());
	for (int i = 0 ; i < kgrid.get_np_red(); ++i)
		reducibleIndices[i] = i;
	std::vector<std::vector<std::complex<float>>> wfctsPerK;
	phonyWave.generate_reducible_grid_wfcts(bandIndices,reducibleIndices,wfctsPerK,npwPerK);

	for ( int ikir = 0 ; ikir < kgrid.get_np_irred(); ++ikir)
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
	boost::filesystem::remove(phonyDir / "WAVECAR");
}

BOOST_AUTO_TEST_CASE( Phony_VASP_Wfct_interpolation )
{
	//We start by reading parameters of other files in this directory for lattice matrices and kpoints
	test::fixtures::MockStartup ms;
	boost::filesystem::path phonyDir = ms.get_data_for_testing_dir() / "phony" / "vasp_sym";

	std::string input = std::string()+
			"root_dir = "+phonyDir.string()+"\n";
	elephon::IOMethods::InputOptions opts;
	ms.simulate_elephon_input(
			(phonyDir / "infile").string(),
			input,
			opts);

	//introduce the VASP data loader
	std::shared_ptr<elephon::IOMethods::VASPInterface> loader =
			std::make_shared< elephon::IOMethods::VASPInterface >(opts);

	elephon::LatticeStructure::LatticeModule lattice;
	loader->read_lattice_structure(phonyDir.string(), lattice);
	elephon::LatticeStructure::Symmetry sym;
	sym.set_reciprocal_space_sym();
	elephon::LatticeStructure::RegularSymmetricGrid kgrid;
	kgrid.initialize(
			std::vector<int>{6, 6, 6},
			1e-6,
			std::vector<double>{0.0, 0.0, 0.0},
			sym,
			lattice);

	const int nBnd = 2;
	const double eCut = 11;
	std::vector<double> bandData(kgrid.get_np_irred()*nBnd, 0.0);
	elephon::ElectronicStructure::ElectronicBands bands;
	bands.initialize(nBnd, 0.0, bandData, kgrid);

	elephon::IOMethods::ReadVASPWaveFunction reader;
	std::vector<int> mFFT;
	reader.compute_fourier_max(eCut, lattice, mFFT);
	std::vector<std::vector<std::complex<float>>> wavefunctions(kgrid.get_np_irred());
	std::vector<std::vector<int>> fftMap;
	std::vector<double> kvectors(3*kgrid.get_np_irred());
	for (int ikir = 0 ; ikir < kgrid.get_np_irred(); ++ikir)
	{
		auto k = kgrid.get_vector_direct(kgrid.get_maps_irreducible_to_reducible()[ikir][kgrid.get_symmetry().get_identity_index()]);
		kvectors[ikir*3+0] = k[0];
		kvectors[ikir*3+1] = k[1];
		kvectors[ikir*3+2] = k[2];
	}
	reader.compute_fourier_map(
			kvectors,
			fftMap,
			kgrid.get_grid_prec(),
			1,
			mFFT,
			eCut,
			lattice);

	std::vector<int> npwK(kgrid.get_np_irred());
	for (int ikir = 0 ; ikir < kgrid.get_np_irred(); ++ikir)
	{
		npwK[ikir] = fftMap[ikir].size()/3;
		wavefunctions[ikir].resize(npwK[ikir]*nBnd);
		for ( int ipw = 0 ; ipw < npwK[ikir]; ++ipw)
		{
			int ibnd = 0;
			wavefunctions[ikir][ibnd*npwK[ikir] + ipw] =
					std::complex<float>(kvectors[ikir*3+0], kvectors[ikir*3+2]);
			ibnd = 1;
			wavefunctions[ikir][ibnd*npwK[ikir] + ipw] = std::complex<float>(kvectors[ikir*3+1], ibnd);
		}
	}

	elephon::IOMethods::write_VASP_wavefunctions(
			(phonyDir / "WAVECAR" ).string(),
			eCut,
			lattice,
			bands,
			wavefunctions,
			fftMap);

	std::vector<double> k{0.0, 0.0, 0.0,
						  1.0/3.0, 0.275, 0.0,
						  0.0, 0.1, 0.025,};
	std::vector<int> bandsList{0, 1};

	elephon::ElectronicStructure::Wavefunctions wfcts;
	wfcts.initialize(
			kgrid,
			loader);

	std::vector<std::vector<std::complex<float>>> wfctsArbK;
	std::vector<std::vector<int>> fftMapArbK;
	wfcts.generate_wfcts_at_arbitray_kp(
	                k,
	                bandsList,
	                wfctsArbK,
					fftMapArbK);

	BOOST_REQUIRE_EQUAL(int(fftMapArbK.size()), 3);

	// k = 0 0 0
	int ipw = 0; // This must be a plane wave which is not part of the border such that a zero mixes in ...
	int npw = fftMapArbK[0].size()/3;
	int ibnd = 0;
	BOOST_CHECK_SMALL( std::abs(wfctsArbK[0][ibnd*npw + ipw]
								- std::complex<float>(0.0, 0.0)), float(1e-5) );
	ibnd = 1;
	BOOST_CHECK_SMALL( std::abs(wfctsArbK[0][ibnd*npw + ipw]
								- std::complex<float>(0.0, 1.0)), float(1e-5) );

	// k = 0.333... 0.275 0
	npw = fftMapArbK[1].size()/3;
	ibnd = 0;
	BOOST_CHECK_SMALL( std::abs(wfctsArbK[1][ibnd*npw + ipw]
								- std::complex<float>(1.0/3.0, 0.0)), float(1e-5) );
	ibnd = 1;
	BOOST_CHECK_SMALL( std::abs(wfctsArbK[1][ibnd*npw + ipw]
								- std::complex<float>(0.275, 1.0)), float(1e-5) );

	// k = 0.0, 0.1, 0.025
	npw = fftMapArbK[2].size()/3;
	ibnd = 0;
	BOOST_CHECK_SMALL( std::abs(wfctsArbK[2][ibnd*npw + ipw]
								- std::complex<float>(0.0, 0.025)), float(1e-5) );
	ibnd = 1;
	BOOST_CHECK_SMALL( std::abs(wfctsArbK[2][ibnd*npw + ipw]
								- std::complex<float>(0.1, 1.0)), float(1e-5) );

	boost::filesystem::remove(phonyDir / "WAVECAR");
}

BOOST_AUTO_TEST_CASE( MgB2_vasp_wfct_arbitray_kpts )
{
	auto wfcts = get_MgB2_vasp_wfct();
	std::vector<double> k{0.01, 0.01, 0.01};
	auto starOps = wfcts.get_k_grid().get_symmetry().star_operations(k);
	std::vector<double> ks( (starOps.size()+1)*3 );
	std::copy(k.begin(), k.end(), &ks[0]);
	for ( int isym = 0; isym < starOps.size(); ++isym)
	{
		auto kt = k;
		starOps[isym].apply(kt, true);
		std::copy(kt.begin(), kt.end(), &ks[isym*3]);
	}
	std::vector<int> bands{3};
	std::vector<std::vector<std::complex<float>>> wfctsArbK;
	std::vector<std::vector<int>> fftMap;
	wfcts.generate_wfcts_at_arbitray_kp(
			ks,
			bands,
			wfctsArbK,
			fftMap);

	std::vector<std::complex<float>> wfctsGammaFullGrid;
	elephon::Algorithms::FFTInterface fft;
	std::vector<int> chargeDim = {40,40,48};

	std::vector<double> chggam(chargeDim[0]*chargeDim[1]*chargeDim[2], 0.0);
	for ( int isym = 0; isym < (starOps.size()+1); ++isym)
	{

		fft.fft_sparse_data(
				fftMap[isym],
				wfcts.get_max_fft_dims(),
				wfctsArbK[isym],
				bands.size(),
				-1,
				wfctsGammaFullGrid,
				chargeDim,
				false,
				starOps.size()+1);

		for ( int i = 0; i < chggam.size(); ++i )
			chggam[i] += std::pow(std::real(wfctsGammaFullGrid[i]),2)+std::pow(std::imag(wfctsGammaFullGrid[i]),2);
	}

	write_MgB2_vasp_chg(chargeDim, chggam);
}

BOOST_AUTO_TEST_CASE( MgB2_vasp_wfct_Gamma )
{
	auto wfcts = get_MgB2_vasp_wfct();

	std::vector<std::vector<std::complex<float>>> wfctsGamma;
	std::vector<int> npwPerK;
	wfcts.generate_reducible_grid_wfcts(
			std::vector<int>{3},
			std::vector<int>{0},
			wfctsGamma,
			npwPerK);

	std::vector<std::vector<int>> fftMap;
	wfcts.compute_Fourier_maps(
			std::vector<double>{0.0, 0.0, 0.0},
			fftMap);

	BOOST_REQUIRE_EQUAL(npwPerK[0], fftMap[0].size()/3);

	std::vector<int> chargeDim = {40,40,48};

	std::vector<std::complex<float>> wfctsGammaFullGrid;
	elephon::Algorithms::FFTInterface fft;

	fft.fft_sparse_data(
			fftMap[0],
			wfcts.get_max_fft_dims(),
			wfctsGamma[0],
			1,
			-1,
			wfctsGammaFullGrid,
			chargeDim,
			false,
			1);

	double norm2_1 = 0;
	for ( auto a : wfctsGamma[0])
		norm2_1 += std::pow(std::abs(a),2);

	double norm2_2 = 0;
	for ( auto a :wfctsGammaFullGrid)
		norm2_2 += std::pow(std::abs(a),2)/wfctsGammaFullGrid.size();

	BOOST_CHECK_SMALL(norm2_1 - norm2_2, 1e-5);

	std::vector<double> chggam(chargeDim[0]*chargeDim[1]*chargeDim[2]);
	for ( int i = 0; i < chggam.size(); ++i )
		chggam[i] = std::pow(std::real(wfctsGammaFullGrid[i]),2)+std::pow(std::imag(wfctsGammaFullGrid[i]),2);

	write_MgB2_vasp_chg(chargeDim, chggam);
}

BOOST_AUTO_TEST_CASE( MgB2_vasp_wfct_star_k )
{
	// Compute the star of a k point and sum up the charge. It must be symmetric
	auto wfcts = get_MgB2_vasp_wfct();
	elephon::Algorithms::FFTInterface fft;

	//We investigate band 3 + 4 because they are degenerate at the 1/3 1/3 x points
	std::vector<int> chargeDim = {40,40,48};
	auto kg = wfcts.get_k_grid();
	auto rsSym = kg.get_symmetry();
	rsSym.set_reciprocal_space_sym(false);
	std::vector<std::complex<float>> wfctsRSStarFullGrid;
	for ( int ikir = 0 ; ikir < kg.get_np_irred() ; ++ikir)
	{
		std::vector<double> chgStar(chargeDim[0]*chargeDim[1]*chargeDim[2], 0.0);
		std::vector<int> kindicesRed;
		for ( int is = 0; is < kg.get_maps_sym_irred_to_reducible()[ikir].size(); ++is)
			kindicesRed.push_back( kg.get_maps_irreducible_to_reducible()[ikir][is]);
		std::vector<std::vector<std::complex<float>>> wfctsStar;
		std::vector<int> npwPerK;
		wfcts.generate_reducible_grid_wfcts(
				std::vector<int>{3, 4},
				kindicesRed,
				wfctsStar,
				npwPerK);

		for ( int is = 0; is < kg.get_maps_sym_irred_to_reducible()[ikir].size(); ++is)
		{
			std::vector<std::vector<int>> fftMap;
			auto kv = kg.get_vector_direct(kindicesRed[is]);
			wfcts.compute_Fourier_maps(kv, fftMap);

			BOOST_REQUIRE_EQUAL(npwPerK[is], fftMap[0].size()/3);

			fft.fft_sparse_data(
					fftMap[0],
					wfcts.get_max_fft_dims(),
					wfctsStar[is],
					2,
					-1,
					wfctsRSStarFullGrid,
					chargeDim,
					false,
					1);

			for ( int ib = 0; ib < 2; ++ib)
				for ( int i = 0; i < chgStar.size(); ++i )
				{
					chgStar[i] += std::pow(std::real(wfctsRSStarFullGrid[ib*chgStar.size()+i]),2)
								 + std::pow(std::imag(wfctsRSStarFullGrid[ib*chgStar.size()+i]),2);
					assert( chgStar[i] == chgStar[i] );
				}
		}
		auto r = comp_symm_distortion(chgStar, chargeDim, rsSym);
		BOOST_CHECK_SMALL(r, 1e-5);
	}
}

BOOST_AUTO_TEST_CASE( Al_vasp_wfct_star_k )
{
	// Compute the star of a k point and sum up the charge. It must be symmetric
	auto wfcts = load_Al_fcc_vasp_wfcts();
	elephon::Algorithms::FFTInterface fft;

	//We investigate band 3 + 4 because they are degenerate at the 1/3 1/3 x points
	std::vector<int> chargeDim = {32,32,32};
	auto kg = wfcts.get_k_grid();
	auto rsSym = kg.get_symmetry();
	rsSym.set_reciprocal_space_sym(false);
	std::vector<std::complex<float>> wfctsRSStarFullGrid;

	//We don't test all stars, we test only 10 random samples
	srand (time(NULL));
	std::set<int> samples;
	while( samples.size() < 10 )
	{
		samples.insert( rand() % kg.get_np_irred() );
	}

	for ( int ikir : samples )
	{
		std::cout << "Testing ikir = " << ikir << std::endl;
		std::vector<double> chgStar(chargeDim[0]*chargeDim[1]*chargeDim[2], 0.0);
		std::vector<int> kindicesRed;
		for ( int is = 0; is < kg.get_maps_sym_irred_to_reducible()[ikir].size(); ++is)
			kindicesRed.push_back( kg.get_maps_irreducible_to_reducible()[ikir][is]);
		std::vector<std::vector<std::complex<float>>> wfctsStar;
		std::vector<int> npwPerK;
		wfcts.generate_reducible_grid_wfcts(
				std::vector<int>{0},
				kindicesRed,
				wfctsStar,
				npwPerK);

		for ( int is = 0; is < kg.get_maps_sym_irred_to_reducible()[ikir].size(); ++is)
		{
			std::vector<std::vector<int>> fftMap;
			auto kv = kg.get_vector_direct(kindicesRed[is]);
			wfcts.compute_Fourier_maps(kv, fftMap);

			BOOST_REQUIRE_EQUAL(npwPerK[is], fftMap[0].size()/3);

			fft.fft_sparse_data(
					fftMap[0],
					wfcts.get_max_fft_dims(),
					wfctsStar[is],
					1,
					-1,
					wfctsRSStarFullGrid,
					chargeDim,
					false,
					1);

			for ( int ib = 0; ib < 1; ++ib)
				for ( int i = 0; i < chgStar.size(); ++i )
				{
					chgStar[i] += std::pow(std::real(wfctsRSStarFullGrid[ib*chgStar.size()+i]),2)
								 + std::pow(std::imag(wfctsRSStarFullGrid[ib*chgStar.size()+i]),2);
					assert( chgStar[i] == chgStar[i] );
				}
		}
		auto r = comp_symm_distortion(chgStar, chargeDim, rsSym);
		BOOST_CHECK_SMALL(r, 1e-5);
	}
}

