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
	auto bands = resourceHandler->get_dense_electronic_bands_obj();
	auto ph = resourceHandler->get_phonon_obj();
	int nModes = ph->get_num_modes();

	// set the frequency grid
	auto phGrid = resourceHandler->get_phonon_grid_obj();
	auto frequencies =  phGrid->setup_frequency_grid(
			resourceHandler->get_optns().get_phrange(),
			resourceHandler->get_optns().get_phnpts());
	freqMin_ = frequencies.front();
	freqMax_ = frequencies.back();
	freqNPts_ = resourceHandler->get_optns().get_phnpts();
	a2F_.assign(freqNPts_ , 0.0);

	auto dvscf = resourceHandler->get_displacement_potential_obj();
	auto interpolationKMesh = resourceHandler->get_interpol_reci_mesh_obj();

	Algorithms::TrilinearInterpolation trilin(*interpolationKMesh);

	// obtain a Fermi surface (set of constant energy surfaces) as a list of k point and weights
	int nSamples = resourceHandler->get_optns().get_numFS();
	auto equalEnergySurfaces = resourceHandler->get_optns().get_ea2f();

	std::vector<double> DOS;
	bands->compute_DOS_tetra(
			resourceHandler->get_tetrahedra_grid(),
			equalEnergySurfaces,
			DOS);

	// load the isosurface
	auto tetraIso = resourceHandler->get_tetrahedra_isosurface();
	const float eVToTHz = 241.79893;

	for ( int iE = 0 ; iE < equalEnergySurfaces.size() ; ++iE )
	{
		auto e = equalEnergySurfaces[iE];
		std::cout << "Calculating a2F(w) at energy e=" << e << "eV relative to the Fermi level\n";
		std::cout << "\tComputing electron-phonon coupling for\n";

		for (int ibnd = 0 ; ibnd < tetraIso->get_nBnd(); ++ibnd )
		{
			std::vector<double> isoK, isoWeights;
			tetraIso->get_irreducible_iso_vector_integration_weights(iE, ibnd, isoK, isoWeights);
			if ( isoK.empty() )
				continue;

			for (int ibndP = 0 ; ibndP < tetraIso->get_nBnd(); ++ibndP )
			{
				std::vector<double> isoKP, isoWeightsP;
				tetraIso->get_reducible_iso_vector_integration_weights(iE, ibndP, isoKP, isoWeightsP);
				if ( isoKP.empty() )
					continue;

				std::cout << "\tbands ("<< ibnd << ") ==> ("  << ibndP << ")" << std::endl;

				// compute the electron phonon matrix elements between these k points
				ElectronPhononCoupling gkkp;
				gkkp.generate_gkkp_and_phonon(isoK, isoKP, {ibnd}, {ibndP}, ph, dvscf, wfcts);

				// integrate the function on each constant energy surfaces.
				for ( int ikf1 = 0 ; ikf1 < isoWeights.size() ; ++ikf1)
					for ( int ikf2 = 0 ; ikf2 < isoWeightsP.size() ; ++ikf2)
					{
						std::vector<std::complex<float>>::iterator gkkpitB, gkkpitE;
						std::vector<float>::iterator phitB, phitE;
						gkkp.get_local_matrix_range(ikf1, ikf2, gkkpitB, gkkpitE, phitB, phitE);
						assert( (std::distance(gkkpitB, gkkpitE) == nModes)
								&& (std::distance(phitB, phitE) == nModes));

						auto phononFreqIt = phitB;
						std::vector<float> gkkpModSqr(nModes);
						for ( int inu = 0; gkkpitB != gkkpitE; ++gkkpitB, ++inu, ++phononFreqIt )
							gkkpModSqr[inu] = isoWeights[ikf1]*isoWeightsP[ikf2]*std::real((*gkkpitB)*std::conj(*gkkpitB))
												/ 2.0f /(*phononFreqIt) / DOS[iE] * eVToTHz;

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

} /* namespace PhononStructure */
} /* namespace elephon */
