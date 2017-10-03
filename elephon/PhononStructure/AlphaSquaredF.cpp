/*	This file AlphaSquaredF.cpp is part of elephon.
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
 *  Created on: Sep 26, 2017
 *      Author: A. Linscheid
 */

#include "PhononStructure/AlphaSquaredF.h"
#include "PhononStructure/ElectronPhononCoupling.h"
#include "ElectronicStructure/FermiSurface.h"
#include "Algorithms/FFTInterface.h"
#include "Algorithms/TrilinearInterpolation.h"
#include "Algorithms/LocalDerivatives.h"
#include <boost/filesystem.hpp>
#include <fstream>
#include <iostream>
#include <chrono>
#include <ctime>

namespace elephon
{
namespace PhononStructure
{

void
AlphaSquaredF::compute_a2F( std::shared_ptr<IOMethods::ElectronicStructureCodeInterface> dataLoader )
{
	std::string root_dir = dataLoader->get_optns().get_elphd();
	boost::filesystem::path calcRoot( root_dir );
	std::string electrons_dir = dataLoader->get_optns().get_eld();

	// load the wavefunctions
	ElectronicStructure::Wavefunctions wfcts;
	wfcts.initialize( electrons_dir , dataLoader);

	// load the unit cell data
	LatticeStructure::UnitCell unitCell;
	dataLoader->read_unit_cell(electrons_dir, dataLoader->get_optns().get_gPrec(), unitCell);

	// load the band structure
	ElectronicStructure::ElectronicBands bands;
	dataLoader->read_band_structure(electrons_dir, bands);

	// rebuild the supercell.
	auto scd = dataLoader->get_optns().get_scell();
	auto supercell = unitCell.build_supercell(scd[0], scd[1], scd[2]);

	// rebuild the irreducible displacements
	std::vector<elephon::LatticeStructure::AtomDisplacement> irreducibleDispl;
	unitCell.generate_displacements(
			dataLoader->get_optns().get_magdispl(),
			dataLoader->get_optns().get_symDispl(),
			irreducibleDispl);

	// read the forces and build the matrix of force constants
	std::vector<std::vector<double>> forces(irreducibleDispl.size());
	for ( int irrep = 0 ; irrep < irreducibleDispl.size(); ++irrep)
	{
		auto irrepD = calcRoot / (std::string("displ_")+std::to_string(irrep));
		if ( ! boost::filesystem::exists(irrepD) )
			throw std::runtime_error(std::string("Error: Directory ")+irrepD.string()+" not present."
					" Must contain a converged VASP run");
		std::vector<double> f;
		dataLoader->read_forces(irrepD.string(), f);
		forces[irrep]= std::move(f);
	}
	auto irrepPast = calcRoot / (std::string("displ_")+std::to_string(irreducibleDispl.size()));
	if ( boost::filesystem::exists(irrepPast) )
		throw std::runtime_error(std::string("Directory ")+irrepPast.string()+" is present but this run has only"
				+std::to_string(int(irreducibleDispl.size())-1)+" irreducible displacements.\n"
						"For safety reasons, the code is stopping here. Please clean up first.");
	ForceConstantMatrix fc;
	fc.build(unitCell, supercell, irreducibleDispl, forces );

	// compute the phonon spectrum
	Phonon ph;
	std::vector<double> masses;
	masses.reserve(unitCell.get_atoms_list().size());
	for ( auto const & a : unitCell.get_atoms_list() )
		masses.push_back(a.get_mass());
	ph.initialize(std::move(fc), masses);
	int nModes = ph.get_num_modes();

	this->set_freq_grid( dataLoader->get_optns().get_phrange(), dataLoader->get_optns().get_phnpts());

	// compute the displacement potential
	int nIrdDispl = int(irreducibleDispl.size());
	std::vector<std::vector<double>> displPot( nIrdDispl );
	std::vector<double> thisDisplPot;
	std::vector<int> dim;
	// Here, we read the potential from the vasp output
	for ( int idispl = 0 ; idispl < nIrdDispl; ++idispl )
	{
		dataLoader->read_electronic_potential(
				( calcRoot / ("displ_"+std::to_string(idispl)) ).string(),
				dim,
				thisDisplPot);

		displPot[idispl] = std::move(thisDisplPot);
	}
	LatticeStructure::RegularBareGrid rsGridSC;
	rsGridSC.initialize( dim, false, dataLoader->get_optns().get_gPrec(), {0.0, 0.0, 0.0}, supercell.get_lattice());

	//Read the normal periodic potential
	dataLoader->read_electronic_potential(
			dataLoader->get_optns().get_root_dir(),
			dim,
			thisDisplPot);
	elephon::LatticeStructure::RegularBareGrid rsGridUC;
	rsGridUC.initialize( dim, false, dataLoader->get_optns().get_gPrec(), {0.0, 0.0, 0.0}, unitCell.get_lattice());
	DisplacementPotential dvscf;
	dvscf.build(
			unitCell,
			supercell,
			irreducibleDispl,
			std::move(rsGridUC),
			std::move(rsGridSC),
			std::move(thisDisplPot),
			std::move(displPot));

	LatticeStructure::RegularBareGrid interpolationMesh;
	interpolationMesh.initialize(
			bands.get_grid().interpret_fft_dim_input(dataLoader->get_optns().get_fftd()),
			true,
			bands.get_grid().get_grid_prec(),
			dataLoader->get_optns().get_ffts(),
			bands.get_grid().get_lattice());

	Algorithms::TrilinearInterpolation trilin(interpolationMesh);

	// obtain a Fermi surface (set of constant energy surfaces) as a list of k point and weights
	int nSamples = dataLoader->get_optns().get_numFS();
	auto equalEnergySurfaces = dataLoader->get_optns().get_ea2f();
	for ( auto e : equalEnergySurfaces)
	{
		std::cout << "Calculating a2F(w) at energy " << e << "eV"<< std::endl;

		auto bndsCrossing = bands.get_bands_crossing_energy_lvls({e});
		std::vector<double> reducibleData;
		bands.generate_interpolated_reducible_data(
				bndsCrossing,
				interpolationMesh,
				reducibleData);

		ElectronicStructure::FermiSurface fs;
		fs.triangulate_surface(
				interpolationMesh,
				bndsCrossing.size(),
				reducibleData,
				nSamples,
				e);

		std::cout << "Found " << fs.get_Fermi_vectors().size()/3 << " Fermi vectors" <<std::endl;

		for ( int ib1 = 0 ; ib1 < bndsCrossing.size(); ++ib1 )
		{
			std::vector<double> kf1, w1;
			fs.obtain_irreducible_Fermi_vectors_for_band(ib1, bands.get_grid().get_symmetry(), kf1, w1);

			// compute the Fermi velocities at vectors kf1 and build the surface integral weight
			std::vector<double> FermiVelocitiesKf1;
			std::vector<int> requiredGridIndices;
			std::vector<double> gradDataAtRequestedIndices;
			trilin.data_query( kf1, requiredGridIndices );
			auto interpolated_band_data_loader = [&reducibleData, &bndsCrossing, &ib1] (int ikr, int ib){
				// define the data loader that fetches the correct band for given 'ib1'
				// note that compute_derivatives_sqr_polynom will call with ib = 0, since we only want this
				// particular band ib1. However, in terms of bndsCrossing, this can be another band.
				assert( ikr*bndsCrossing.size() + ib1 < reducibleData.size());
				return reducibleData[ikr*bndsCrossing.size() + ib1];
				};
			Algorithms::localDerivatives::compute_derivatives_sqr_polynom<double>(
					1,
					requiredGridIndices,
					&gradDataAtRequestedIndices,
					nullptr,
					interpolationMesh,
					interpolated_band_data_loader);
			assert(gradDataAtRequestedIndices.size() == 3*requiredGridIndices.size());
			trilin.interpolate(3, gradDataAtRequestedIndices, FermiVelocitiesKf1);
			for ( int ikf = 0 ; ikf < w1.size(); ++ikf)
			{
				double modGradE = std::sqrt(std::pow(FermiVelocitiesKf1[ikf*3+0],2)
										+std::pow(FermiVelocitiesKf1[ikf*3+1],2)
										+std::pow(FermiVelocitiesKf1[ikf*3+2],2));
				if ( modGradE < 1e-6) //cutoff
					modGradE = 1e-6;
				w1[ikf] = w1[ikf]/modGradE*interpolationMesh.get_lattice().get_volume()/std::pow(2*M_PI,3);
			}

			for ( int ib2 = 0 ; ib2 < bndsCrossing.size(); ++ib2 )
			{
				std::vector<double> kf2 = fs.get_Fermi_vectors_for_band(ib2);
				std::vector<double> w2 = fs.get_Fermi_weights_for_band(ib2);

				// compute the Fermi velocities at vectors kf2
				std::vector<double> FermiVelocitiesKf2;
				trilin.data_query( kf2, requiredGridIndices );
				auto interpolated_band_data_loader_kp = [&reducibleData, &bndsCrossing, &ib2] (int ikr, int ib){
					// define the data loader that fetches the correct band for given 'ib1'
					// note that compute_derivatives_sqr_polynom will call with ib = 0, since we only want this
					// particular band ib1. However, in terms of bndsCrossing, this can be another band.
					assert( ikr*bndsCrossing.size() + ib2 < reducibleData.size());
					return reducibleData[ikr*bndsCrossing.size() + ib2];
					};
				Algorithms::localDerivatives::compute_derivatives_sqr_polynom<double>(
						1,
						requiredGridIndices,
						&gradDataAtRequestedIndices,
						nullptr,
						interpolationMesh,
						interpolated_band_data_loader_kp);
				assert(gradDataAtRequestedIndices.size() == 3*requiredGridIndices.size());
				trilin.interpolate(3, gradDataAtRequestedIndices, FermiVelocitiesKf2);
				for ( int ikf = 0 ; ikf < w2.size(); ++ikf)
				{
					double modGradE = std::sqrt(std::pow(FermiVelocitiesKf2[ikf*3+0],2)
											+std::pow(FermiVelocitiesKf2[ikf*3+1],2)
											+std::pow(FermiVelocitiesKf2[ikf*3+2],2));
					if ( modGradE < 1e-6) //cutoff
						modGradE = 1e-6;
					w2[ikf] = w2[ikf]/modGradE*interpolationMesh.get_lattice().get_volume()/std::pow(2*M_PI,3);
				}

				// compute the electron phonon matrix elements between these k points
				ElectronPhononCoupling gkkp;
				gkkp.generate_gkkp_and_phonon(kf1, kf2, {ib1}, {ib2}, ph, dvscf, wfcts);

				// integrate the function on each constant energy surfaces.
				for ( int ikf1 = 0 ; ikf1 < w1.size() ; ++ikf1)
					for ( int ikf2 = 0 ; ikf2 < w2.size() ; ++ikf2)
					{
						std::vector<std::complex<float>>::iterator gkkpitB, gkkpitE;
						std::vector<float>::iterator phitB, phitE;
						gkkp.get_local_matrix_range(ikf1, ikf2, gkkpitB, gkkpitE, phitB, phitE);
						assert( (std::distance(gkkpitB, gkkpitE) == nModes)
								&& (std::distance(phitB, phitE) == nModes));

						std::vector<float> gkkpModSqr(nModes);
						for ( int inu = 0; gkkpitB != gkkpitE; ++gkkpitB, ++inu )
							gkkpModSqr[inu] = w1[ikf1]*w2[ikf2]*std::real((*gkkpitB)*std::conj(*gkkpitB));

						this->map_freq_grid_slot(gkkpModSqr.begin(), phitB, phitE);
					}
			}
		}
	}

	auto DeltaOmega = (freqMax_ - freqMin_)/freqNPts_;
	for ( auto &w : a2F_)
		w /= DeltaOmega;
}

void
AlphaSquaredF::map_freq_grid_slot(
		std::vector<float>::const_iterator begData,
		std::vector<float>::const_iterator begFreq,
		std::vector<float>::const_iterator endFreq)
{
	assert((a2F_.size() == freqNPts_) && (freqNPts_ > 0));
	auto itD = begData;
	for (auto it = begFreq ; it != endFreq; ++it  )
	{
		double omega = *it;
		int iomega = std::floor((omega-freqMin_)*freqNPts_/(freqMax_ - freqMin_));
		if ( (iomega < 0) or (iomega >= freqNPts_) )
			continue;
		a2F_[iomega] += *itD;
	}
}

void
AlphaSquaredF::write_a2F_file(std::string const & filename) const
{
	std::ofstream file( filename.c_str() );
	if ( not file.good() )
		throw std::runtime_error(std::string("Problem opening file ")+filename+" for writing the a2F data");

	auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
	file << "# a2F(w) with frequency in units of THz. Date is " << std::ctime(&now) << std::endl;

	for ( int iw = 0 ; iw < freqNPts_ ; ++iw)
	{
		double w = freqMin_ + (iw + 0.5)*(freqMax_ - freqMin_)/freqNPts_;
		file << w << '\t' << a2F_[iw] << '\n';
	}
}

void
AlphaSquaredF::set_freq_grid(std::vector<double> const & freqRange, int npts)
{
	if ( freqRange.empty() )
	{
		// determine the phonon range automatically
		throw std::runtime_error("Automatic determination of the phonon range not implemented");
	}
	else if (freqRange.size() == 1)
	{
		freqMin_ = 0;
		freqMax_ = freqRange.front();
	}
	else if (freqRange.size() == 2)
	{
		freqMin_ = freqRange[0];
		freqMax_ = freqRange[1];
	}
	freqNPts_ = npts;
	a2F_.resize(freqNPts_);
}

} /* namespace PhononStructure */
} /* namespace elephon */
