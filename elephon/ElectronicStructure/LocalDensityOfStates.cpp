/*	This file LocalDensityOfStates.cpp is part of elephon.
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
 *  Created on: Jul 24, 2017
 *      Author: A. Linscheid
 */

#include "ElectronicStructure/LocalDensityOfStates.h"
#include "ElectronicStructure/FermiSurface.h"
#include "Algorithms/TrilinearInterpolation.h"
#include "LatticeStructure/UnitCell.h"
#include "Algorithms/FFTInterface.h"
#include "IOMethods/WriteVASPRealSpaceData.h"

namespace elephon
{
namespace ElectronicStructure
{

void
LocalDensityOfStates::compute_ldos(
		std::vector<double> energies,
		Wavefunctions const& wfcts,
		LatticeStructure::UnitCell unitcell,
		int nkpointsPerSurface,
		std::vector<int> realSpaceRes,
		ElectronicBands const & bands)
{
	uc_ = std::move(unitcell);
	isoEnergies_ = std::move(energies);
	rsDims_ = std::move(realSpaceRes);
	auto reqBandIds = bands.get_bands_crossing_energy_lvls(energies);
	std::vector<double> regularData;
	bands.generate_reducible_grid_bands(reqBandIds, regularData);

	GradientFFTReciprocalGrid gradE;
	gradE.compute_gradient(
			bands.get_grid().get_grid_dim(),
			uc_.get_lattice(),
			reqBandIds.size(),
			regularData);

	Algorithms::TrilinearInterpolation triLin( wfcts.get_k_grid().view_bare_grid() );
	Algorithms::FFTInterface fft;

	std::vector< std::vector< std::complex<float> > > wfctsFs;
	std::vector< std::complex<float> > wfctsRealSpace;
	std::vector< std::vector< int > > fftMapsWfctsFs;
	std::vector<double> gradDataAtRequestedIndices;
	std::vector<int> reqestedIndices;
	std::vector<double> FermiVelocities;
	std::vector< std::complex<float> > wfctsOneBand;
	int nrs = rsDims_[0]*rsDims_[1]*rsDims_[2];
	ldos_.assign( isoEnergies_.size()*nrs , 0.0 );
	for ( int ie = 0 ; ie < isoEnergies_.size(); ++ie )
	{
		double e = isoEnergies_[ie];
		ElectronicStructure::FermiSurface fs;
		fs.triangulate_surface(
				bands.get_grid().get_grid_dim(),
				uc_.get_lattice(),
				reqBandIds.size(),
				regularData,
				nkpointsPerSurface,
				e);

		wfcts.generate_wfcts_at_arbitray_kp(
				fs.get_Fermi_vectors(),
				reqBandIds,
				wfctsFs,
				fftMapsWfctsFs);

		for ( int ibRel = 0 ; ibRel < reqBandIds.size(); ++ibRel )
		{
			int bandOffset = fs.get_band_offset(ibRel);
			auto const & kfv = fs.get_Fermi_vectors_for_band(ibRel);
			auto const & kfw = fs.get_Fermi_weights_for_band(ibRel);
			if ( kfv.size() == 0 )
				continue;

			//Compute the Fermi velocities
			triLin.data_query( kfv, reqestedIndices );
			gradE.copy_data( reqestedIndices, std::vector<int>{ibRel}, gradDataAtRequestedIndices );
			triLin.interpolate(3,gradDataAtRequestedIndices,FermiVelocities);

			for (int ikf = 0 ; ikf < kfv.size()/3; ++ikf)
			{
				int npw = fftMapsWfctsFs[bandOffset+ikf].size()/3;
				wfctsOneBand.resize(npw);
				std::copy(&wfctsFs[bandOffset+ikf][ibRel*npw], &wfctsFs[bandOffset+ikf][ibRel*npw] + npw,
						wfctsOneBand.begin());

				//Compute the wavefunction part |psi(r)|^2
				fft.fft_sparse_data(
						fftMapsWfctsFs[bandOffset+ikf],
						wfcts.get_max_fft_dims(),
						wfctsOneBand,
						1,
						-1,
						wfctsRealSpace,
						rsDims_,
						false, //we want C data order
						nkpointsPerSurface);

				double modGradE = std::sqrt(std::pow(FermiVelocities[ikf*3+0],2)
										+std::pow(FermiVelocities[ikf*3+1],2)
										+std::pow(FermiVelocities[ikf*3+2],2));
				if ( modGradE < 1e-6) //cutoff
					modGradE = 1e-6;
				double contrib = kfw[ikf]/modGradE*unitcell.get_lattice().get_volume()/std::pow(2*M_PI,3);
				for ( int ir = 0 ; ir < nrs; ++ir )
				{
					ldos_[ie*nrs+ir] += contrib*(std::pow(std::real(wfctsRealSpace[ir]), 2)
										       + std::pow(std::imag(wfctsRealSpace[ir]), 2));
					assert( ! std::isnan(ldos_[ie*nrs+ir]) );
				}
			}
		}
	}
}

void
LocalDensityOfStates::compute_ldos(
		std::vector<double> const & energies,
		std::shared_ptr<IOMethods::ElectronicStructureCodeInterface> loader )
{
	ElectronicBands bands;
	loader->read_band_structure(
			loader->get_optns().get_root_dir(),
			bands);

	LatticeStructure::RegularSymmetricGrid kgrid;
	LatticeStructure::LatticeModule lattice;
	std::vector<LatticeStructure::Atom> atoms;
	LatticeStructure::Symmetry sym;
	loader->read_cell_paramters(
			loader->get_optns().get_root_dir(),
			loader->get_optns().get_gPrec(),
			kgrid,
			lattice,
			atoms,
			sym);
	LatticeStructure::UnitCell uc;
	uc.initialize(atoms, lattice, sym);

	Wavefunctions wfcts;
	wfcts.initialize(
			loader->get_optns().get_root_dir(),
			loader);

	//Charge grid in general has to be 2x the FFT grid of the wfcts
	auto rsgrid = loader->get_max_fft_dims();
	for ( auto &xi : rsgrid)
		xi *= 2;

	this->compute_ldos(
			energies,
			wfcts,
			uc,
			loader->get_optns().get_numFS(),
			rsgrid,
			bands);
}

void
LocalDensityOfStates::write_file( std::string const & filename, bool binary ) const
{
	IOMethods::WriteVASPRealSpaceData output;
	auto filename_mod = filename;

	int nrs = rsDims_[0]*rsDims_[1]*rsDims_[2];
	for ( int ie = 0 ; ie < isoEnergies_.size(); ++ie )
	{
		if ( isoEnergies_.size() > 1 )
			filename_mod = "e"+std::to_string(ie)+"_"+filename;

		output.write_file(
				filename_mod,
				std::string("LDOS at E=")+std::to_string(isoEnergies_[ie]),
				rsDims_,
				uc_,
				std::vector<double>(&ldos_[ie*nrs],&ldos_[ie*nrs]+nrs),
				false,
				true);
	}
}

} /* namespace ElectronicStructure */
} /* namespace elephon */