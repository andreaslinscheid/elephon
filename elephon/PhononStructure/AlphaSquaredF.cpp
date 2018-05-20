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
#include "ElectronicStructure/TetrahedraIsosurface.h"
#include "Algorithms/FFTInterface.h"
#include "Algorithms/TrilinearInterpolation.h"
#include "Algorithms/LocalDerivatives.h"
#include "Auxillary/UnitConversion.h"
#include "PhononStructure/PhononGrid.h"
#include "Algorithms/SimpsonIntegrator.h"
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
AlphaSquaredF::compute_a2F_grid( std::shared_ptr<IOMethods::ResourceHandler> resourceHandler )
{
	auto wfcts = resourceHandler->get_wfct_obj();
	auto bands = resourceHandler->get_dense_electronic_bands_obj();
	auto ph = resourceHandler->get_phonon_obj();
	int nModes = ph->get_num_modes();

	// set the frequency grid
	auto phGrid = resourceHandler->get_phonon_grid_obj();
	this->setup_internal_freq_grid(
			phGrid,
			resourceHandler->get_optns().get_phrange(),
			resourceHandler->get_optns().get_phnpts());

	auto dvscf = resourceHandler->get_displacement_potential_obj();
	auto interpolationKMesh = resourceHandler->get_interpol_reci_tetra_mesh_obj();

	Algorithms::TrilinearInterpolation trilin(interpolationKMesh);

	// obtain a Fermi surface (set of constant energy surfaces) as a list of k point and weights
	auto equalEnergySurfaces = resourceHandler->get_optns().get_ea2f();

	std::vector<double> DOS;
	bands->compute_DOS_tetra(
			resourceHandler->get_tetrahedra_grid(),
			equalEnergySurfaces,
			DOS);

	// load the isosurface
	auto tetraIso = resourceHandler->get_tetrahedra_isosurface();
	ElectronPhononCoupling gkkp;

	Auxillary::alignedvector::FV gkkpMod2;
	for ( int iE = 0 ; iE < equalEnergySurfaces.size() ; ++iE )
	{
		auto e = equalEnergySurfaces[iE];
		std::cout << "Calculating a2F(w) at energy e=" << e << "eV relative to the Fermi level\n"
				"using matrix elements on the dense electronic grid\n";

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

				const int nKF = isoWeights.size();
				const int nKFP = isoWeightsP.size();

				std::cout << "\tbands ("<< ibnd << ") ==> ("  << ibndP << ")" << std::endl;

				// perform a interpolation of the phonon frequency grid
				std::vector<int> requestPhononGridIndices;
				std::vector<double> qfvec(3*nKF*nKFP);
				for (int ikf = 0 ; ikf < nKF; ++ikf)
					for (int ikfp = 0 ; ikfp < nKFP; ++ikfp)
						for (int i = 0 ; i < 3; ++i)
							qfvec[(ikf*nKFP+ikfp)*3+i] = isoKP[ikfp*3+i]-isoK[ikf*3+i];
				trilin.data_query(qfvec, requestPhononGridIndices);

				std::vector<double> isoPHData;
				std::vector<double> phononGridData;
				phGrid->copy_selected_grid_points(requestPhononGridIndices, phononGridData);
				trilin.interpolate( nModes,
									phononGridData,
									isoPHData);
				const int nqA = qfvec.size()/3;

				// get a list of grid-q points that are necessary.
				// also compute the list of k- and k'- grid points that are the closes to the
				// iso-surface k and k' point (qGrid = k'Grid - kGrid ).
				std::vector<std::pair<int,
									  std::vector<std::pair<int,int>>>> q_to_k_and_kp_map;
				std::vector<std::pair<int,std::vector<int>>> kGridToKList,
															kpGridToKpList;
				this->query_q(
						isoK,
						isoKP,
						kGridToKList,
						kpGridToKpList,
						wfcts->get_k_grid().view_bare_grid(),
						q_to_k_and_kp_map);

				// for each such grid-q point compute gkkp and integrate all k and k' iso-surface functions
				// connected with this grid-q point.
				int nEval = 0;
				for (auto const & gridQ : q_to_k_and_kp_map)
				{
					// gridQ.second contains the indices in the kGridToKList and kpGridToKpList lists,
					// here, we need to convert it to the actual reducible k grid indices for the gkkp
					// matrix element calculation.
					auto kAndKPGridIndesList = gridQ.second;
					for (auto & p : kAndKPGridIndesList)
					{
						p.first = kGridToKList[p.first].first;
						p.second = kpGridToKpList[p.second].first;
					}

					gkkp.generate_gkkp_mod_2_of_q(
							gridQ.first,
							kAndKPGridIndesList,
							{ibnd},
							{ibndP},
							ph,
							dvscf,
							wfcts,
							gkkpMod2);

					for (int ikThisQ = 0 ; ikThisQ < gridQ.second.size(); ++ikThisQ)
					{
						const int ikGridList = gridQ.second[ikThisQ].first;
						const int ikpGridList = gridQ.second[ikThisQ].second;

						assert(((ikGridList >=0) && (ikGridList<kGridToKList.size()))
								&& ((ikpGridList >=0) && (ikpGridList<kpGridToKpList.size())));
						for (int ikf : kGridToKList[ikGridList].second )
							for (int ikfp : kpGridToKpList[ikpGridList].second )
							{
								assert((ikf>=0)&&(ikf<isoWeights.size()));
								assert((ikfp>=0)&&(ikfp<isoWeightsP.size()));

								float w1 = isoWeights[ikf];
								float w2 = isoWeightsP[ikfp];
								nEval++;
								for ( int inu = 0; inu < nModes; ++inu)
								{
									double omega = isoPHData[(ikf*nKFP+ikfp)*nModes+inu];
									int iomega = std::floor((omega-freqMin_)*freqNPts_/(freqMax_ - freqMin_));
									if ( (iomega < 0) or (iomega >= freqNPts_) )
										continue;
									a2F_[iomega] += static_cast<double>(w1*w2*gkkpMod2[ikThisQ*nModes+inu]) / DOS[iE];
								}
							}
					}
				}
				if(nEval!=nqA)
					throw std::logic_error("Incorrect number of evaluations.");
			}
		}
	}

	auto DeltaOmega = (freqMax_ - freqMin_)/freqNPts_;
	for ( auto &w : a2F_)
		w *= Auxillary::units::EV_TO_THZ_CONVERSION_FACTOR / DeltaOmega;
}

void
AlphaSquaredF::setup_internal_freq_grid(
		std::shared_ptr<const PhononStructure::PhononGrid> phgrid,
		std::vector<double> const & phrange,
		int npts)
{
	freqNPts_ = npts;
	frequencies_.resize(npts);
	auto frequencies = phgrid->setup_frequency_grid(phrange, freqNPts_);
	std::copy(frequencies.begin(),frequencies.end(), frequencies_.begin());
	freqMin_ = frequencies_.front();
	freqMax_ = frequencies_.back();
	a2F_.assign(freqNPts_ , 0.0);
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
AlphaSquaredF::load_from_file(std::string const & filename)
{
	std::ifstream file( filename.c_str() );
	if ( not file.good() )
		throw std::runtime_error(std::string("Problem opening file ")+filename+" for reading the a2F data");

	frequencies_.clear();
	a2F_.clear();

	std::string s;
	std::getline(file,s); // comment line
	numBands_ = -1;
	std::vector<double> data;
	while( std::getline(file,s) )
	{
		std::stringstream ss(s);
		std::vector<double> values;
		double element;
		while ( ss >> element ){
			values.push_back(element);
		}
		if ( values.empty())
			continue;

		if ( numBands_ < 0 )
			numBands_ = sqrt(values.size()-1);

		if ( (numBands_*numBands_ != values.size()-1) && (values.size() < 2))
			throw std::runtime_error(std::string("Problem reading file ")+filename+" for reading the a2F data\n"
					"Each non-empty line must have the same number of values, frequency + list of bands data");
		frequencies_.push_back(*values.begin());
		if ( (frequencies_.size() > 1) && (*frequencies_.rbegin() <= *(++frequencies_.rbegin())) )
				throw std::runtime_error(std::string("Problem reading file ")+filename+" for reading the a2F data\n"
						"The frequencies are not in increasing order.");
		data.insert(data.end(), values.begin()+1, values.end());
	}
	frequencies_.shrink_to_fit();
	freqNPts_ = frequencies_.size();
	a2F_.resize(freqNPts_*numBands_*numBands_);

	for (int iband = 0 ; iband < numBands_; ++iband)
		for (int ibandPrime = 0 ; ibandPrime < numBands_; ++ibandPrime)
			for (int ifreq = 0 ; ifreq < freqNPts_; ++ifreq)
				a2F_[ifreq + freqNPts_*(iband + numBands_*ibandPrime)] = data[iband + numBands_*(ibandPrime+numBands_*ifreq)];

	freqMin_ = *frequencies_.begin();
	freqMax_ = *frequencies_.rbegin();
}

void
AlphaSquaredF::query_q(
		std::vector<double> const & kList,
		std::vector<double> const & kpList,
		std::vector<std::pair<int,std::vector<int>>> & kPointToGrid,
		std::vector<std::pair<int,std::vector<int>>> & kpPointToGrid,
		LatticeStructure::RegularBareGrid const & grid,
		std::vector<
					std::pair<int, std::vector<std::pair<int,int>>>> & q_to_k_and_kp_map)
{
	assert(kList.size()%3 == 0);
	assert(kpList.size()%3 == 0);
	const int nK = kList.size()/3;
	const int nKp = kpList.size()/3;

	if (nK*nKp == 0)
		return;

	// NOTE: q = k' - k
	std::vector<int> closestGPK, closestGPKp;
	grid.find_closest_reducible_grid_points(kList, closestGPK);
	grid.find_closest_reducible_grid_points(kpList, closestGPKp);

	std::map<int,std::vector<int>> kPointToGridMap;
	std::map<int,std::vector<int>> kpPointToGridMap;
	for (int ikArbitrayPoint = 0 ; ikArbitrayPoint < closestGPK.size(); ++ikArbitrayPoint )
	{
		auto ret = kPointToGridMap.insert(std::make_pair(closestGPK[ikArbitrayPoint], std::vector<int>()));
		ret.first->second.push_back(ikArbitrayPoint);
	}

	for (int ikArbitrayPoint = 0 ; ikArbitrayPoint < closestGPKp.size(); ++ikArbitrayPoint )
	{
		auto ret = kpPointToGridMap.insert(std::make_pair(closestGPKp[ikArbitrayPoint], std::vector<int>()));
		ret.first->second.push_back(ikArbitrayPoint);
	}

	const int nKGrid = kPointToGridMap.size();
	const int nKpGrid = kpPointToGridMap.size();

	kPointToGrid.clear();
	kPointToGrid.reserve(nKGrid);
	for (auto const & kg : kPointToGridMap)
		kPointToGrid.push_back(kg);

	kpPointToGrid.clear();
	kpPointToGrid.reserve(nKpGrid);
	for (auto const & kg : kpPointToGridMap)
		kpPointToGrid.push_back(kg);

	const std::int64_t nQ =  nKGrid*nKpGrid;

	std::vector<double> allQVectors(nKGrid*nKpGrid*3);
	std::vector<std::pair<int,int>> q_index_to_k_and_kp_map(nKGrid*nKpGrid);
	std::vector<double> k(3), kp(3);
	std::vector<int> xyzBuffer(3);
	for ( int ik = 0 ; ik < nKGrid; ++ik)
		for ( int ikp = 0 ; ikp < nKpGrid; ++ikp)
		{
			const std::int64_t ikkp = ik*nKpGrid+ikp;
			grid.get_vector_direct(kpPointToGrid[ikp].first, xyzBuffer, kp);
			grid.get_vector_direct(kPointToGrid[ik].first, xyzBuffer, k);
			allQVectors[ikkp*3+0] = kp[0]-k[0];
			allQVectors[ikkp*3+1] = kp[1]-k[1];
			allQVectors[ikkp*3+2] = kp[2]-k[2];

			q_index_to_k_and_kp_map[ikkp].first = ik;
			q_index_to_k_and_kp_map[ikkp].second = ikp;
		}

	// map q vectors back to the 1. BZ
	for ( auto &qi : allQVectors )
		qi -= std::floor(qi+0.5);

	std::vector<int> qGridIndices;
	grid.get_list_reducible_lattice_point_indices(allQVectors, qGridIndices);

	// collect all arbitrary q points for grid points to which they are closest.
	// this results in a mapping from a given grid point to the list of q point
	// indices which are closest to it within the grid.
	std::map<int, std::vector<std::int64_t>> qCollector;
	for (std::int64_t iq = 0 ; iq < nQ ;++iq)
	{
		auto ret = qCollector.insert(std::move(std::make_pair(qGridIndices[iq], std::vector<std::int64_t>())));
		ret.first->second.push_back(iq);
	}

#ifndef NDEBUG
	// cross check that all q-indices are there exactly once.
	std::int64_t nElem = 0;
	std::set<std::int64_t> cross_check_set;
	for (auto gp : qCollector )
	{
		cross_check_set.insert(gp.second.begin(), gp.second.end());
		nElem += gp.second.size();
	}
	assert(nElem == nQ);
	assert(cross_check_set.size() == nQ);
	assert(*cross_check_set.rbegin() == (nQ-1));
#endif

	// now use the q index to k and k' index mapping to resolve the q index in the above data structure
	q_to_k_and_kp_map.reserve(qCollector.size());
	for ( auto const & coarseGP : qCollector )
	{
		auto iqAndKKPListPair = std::make_pair(	coarseGP.first,
												std::vector<std::pair<int,int>>());
		iqAndKKPListPair.second.reserve(coarseGP.second.size());
		for (auto iq : coarseGP.second)
			iqAndKKPListPair.second.push_back(q_index_to_k_and_kp_map[iq]);
		q_to_k_and_kp_map.push_back(std::move(iqAndKKPListPair));
	}
}

double
AlphaSquaredF::operator() (int iFrequency, int iBand, int iBandPrime) const
{
	assert((iFrequency >=0)&&(iFrequency<freqNPts_));
	assert((iBand >=0)&&(iBand<numBands_));
	assert((iBandPrime >=0)&&(iBandPrime<numBands_));
	return a2F_[iFrequency + freqNPts_*(iBand+numBands_*iBandPrime)];
}

int
AlphaSquaredF::get_num_bands() const
{
	return numBands_;
}

Auxillary::alignedvector::DV::const_iterator
AlphaSquaredF::beginBand(int iband, int ibandPrime) const
{
	return a2F_.begin()+freqNPts_*(iband+numBands_*ibandPrime);
}

Auxillary::alignedvector::DV::const_iterator
AlphaSquaredF::endBand(int iband, int ibandPrime) const
{
	return a2F_.begin()+freqNPts_*(iband+numBands_*ibandPrime) + freqNPts_;
}

int
AlphaSquaredF::get_num_frequency_samples() const
{
	return freqNPts_;
}

Auxillary::alignedvector::DV const &
AlphaSquaredF::get_frequencies() const
{
	return frequencies_;
}

double AlphaSquaredF::get_max_frequency_range() const
{
	return freqMax_;
}

double
AlphaSquaredF::compute_coupling_lambda() const
{
	Algorithms::SimpsonIntegrator<decltype(frequencies_)> si;
	double lambda = 0.0;
	Auxillary::alignedvector::DV buffer(freqNPts_);
	for (int iband = 0; iband < numBands_; ++iband)
		for (int ibandPrime = 0; ibandPrime < numBands_; ++ibandPrime)
		{
			for (int ifreq = 0; ifreq < freqNPts_; ++ifreq)
				buffer[ifreq] = 2.0 * (*this)(ifreq, iband, ibandPrime)	/ frequencies_[ifreq];
			lambda = std::max(lambda, si.integrate(buffer.begin(), buffer.end(), frequencies_.begin(), frequencies_.end()));
		}
	return lambda;
}

double
AlphaSquaredF::compute_omegaLog(double lambda) const
{
	if ( lambda == 0.0 )
		lambda = this->compute_coupling_lambda ();
	Algorithms::SimpsonIntegrator<decltype(frequencies_)> si;
	double omega_log = 0.0;
	Auxillary::alignedvector::DV buffer(freqNPts_);
	for (int iband = 0; iband < numBands_; ++iband)
		for (int ibandPrime = 0; ibandPrime < numBands_; ++ibandPrime)
		{
			for (int ifreq = 0; ifreq < freqNPts_; ++ifreq)
				buffer[ifreq] = 2.0/lambda *std::log(frequencies_[ifreq]) * (*this)(ifreq, iband, ibandPrime)	/ frequencies_[ifreq];
			omega_log = std::max(omega_log, std::exp(si.integrate(buffer.begin(), buffer.end(), frequencies_.begin(), frequencies_.end())));
		}
	return omega_log;
}

} /* namespace PhononStructure */
} /* namespace elephon */
