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
AlphaSquaredF::compute_a2F( std::shared_ptr<IOMethods::ResourceHandler> resourceHandler )
{
	auto wfcts = resourceHandler->get_wfct_obj();
	auto bands = resourceHandler->get_electronic_bands_obj();
	auto ph = resourceHandler->get_phonon_obj();
	int nModes = ph->get_num_modes();

	this->set_freq_grid( resourceHandler->get_optns().get_phrange(), resourceHandler->get_optns().get_phnpts());

	auto dvscf = resourceHandler->get_displacement_potential_obj();
	auto interpolationKMesh = resourceHandler->get_interpol_reci_mesh_obj();

	Algorithms::TrilinearInterpolation trilin(*interpolationKMesh);

	// obtain a Fermi surface (set of constant energy surfaces) as a list of k point and weights
	int nSamples = resourceHandler->get_optns().get_numFS();
	auto equalEnergySurfaces = resourceHandler->get_optns().get_ea2f();
	for ( auto e : equalEnergySurfaces)
	{
		std::cout << "Calculating a2F(w) at energy " << e << "eV"<< std::endl;

		auto bndsCrossing = bands->get_bands_crossing_energy_lvls({e});
		std::vector<double> reducibleData;
		bands->generate_interpolated_reducible_data(
				bndsCrossing,
				*interpolationKMesh,
				reducibleData);

		ElectronicStructure::FermiSurface fs;
		fs.triangulate_surface(
				*interpolationKMesh,
				bndsCrossing.size(),
				reducibleData,
				nSamples,
				e);

		std::cout << "Found " << fs.get_Fermi_vectors().size()/3 << " Fermi vectors" <<std::endl;

		for ( int ib1 = 0 ; ib1 < bndsCrossing.size(); ++ib1 )
		{
			std::vector<double> kf1, w1;
			fs.obtain_irreducible_Fermi_vectors_for_band(ib1, bands->get_grid().get_symmetry(), kf1, w1);

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
					*interpolationKMesh,
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
				w1[ikf] = w1[ikf]/modGradE*interpolationKMesh->get_lattice().get_volume()/std::pow(2*M_PI,3);
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
						*interpolationKMesh,
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
					w2[ikf] = w2[ikf]/modGradE*interpolationKMesh->get_lattice().get_volume()/std::pow(2*M_PI,3);
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
