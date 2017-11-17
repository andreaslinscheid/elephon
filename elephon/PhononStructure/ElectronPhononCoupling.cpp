/*	This file ElectronPhononCouling.cpp is part of elephon.
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
 *  Created on: Jul 4, 2017
 *      Author: A. Linscheid
 */

#include "PhononStructure/ElectronPhononCoupling.h"
#include "Algorithms/TrilinearInterpolation.h"
#include "Algorithms/FFTInterface.h"
#include "Auxillary/UnitConversion.h"
#include <set>
#include <map>
#include <chrono>
#include <sstream>
#include <omp.h>

namespace elephon
{
namespace PhononStructure
{

void
ElectronPhononCoupling::generate_gkkp_and_phonon(
		std::vector<double> kList,
		std::vector<double> kpList,
		std::vector<int> bandList,
		std::vector<int> bandpList,
		std::shared_ptr<const Phonon> ph,
		std::shared_ptr<const DisplacementPotential> dvscf,
		std::shared_ptr<const ElectronicStructure::Wavefunctions> wfcts)
{
	assert(kList.size()%3 == 0);
	assert(kpList.size()%3 == 0);
	nK_ = kList.size()/3;
	nKp_ = kpList.size()/3;
	if ( (nK_ == 0) or (nKp_ == 0) )
		return;

	bool report_timings = nKp_ > 100;

	// TODO a system-wide meaure of timings would be desirable
	std::map<std::string,std::pair< decltype(std::chrono::system_clock::now()), std::chrono::duration<double> >> durations;
	auto start_m_c = [&durations] (std::string const & tag) {
		auto it = durations.find(tag) ;
		if ( it == durations.end())
			durations[tag] = std::make_pair(std::chrono::system_clock::now(),
											std::chrono::duration<double>(0.0));
		else
			it->second.first = std::chrono::system_clock::now();
	};
	auto stop_m_c = [&durations] (std::string const & tag) {
		auto it = durations.find(tag) ;
		assert ( it != durations.end());
		auto time_now = std::chrono::system_clock::now();
		it->second.second += std::chrono::duration<double>(time_now - it->second.first);
	};

	if ( report_timings )
		start_m_c("wfcts_gen");

	//Wave functions are on a regular grid and in G (reciprocal) space
	std::vector< std::vector<std::complex<float> > > wfcsk, wfcskp;
	std::vector< std::vector<int> > fftMapsK, fftMapsKp;
	wfcts->generate_wfcts_at_arbitray_kp( kList, bandList, wfcsk, fftMapsK);
	wfcts->generate_wfcts_at_arbitray_kp( kpList, bandpList, wfcskp, fftMapsKp );

	// conjugate the left wfct of the scalar product
	for ( auto & wfk : wfcsk )
		for ( auto & cg : wfk )
			cg = std::conj(cg);

	if ( report_timings )
	{
		stop_m_c("wfcts_gen");
		std::cout << "\tGenerating wavefunctions at "<<kList.size()/3+kpList.size()/3
				<< " k points: "<<durations["wfcts_gen"].second.count()<<"s"<< std::endl;
		durations.erase("wfcts_gen");
	}

	std::vector<int> potentialFFTGrid = dvscf->get_real_space_grid().get_grid_dim();

	int nr = potentialFFTGrid[0]*potentialFFTGrid[1]*potentialFFTGrid[2];
	std::vector<double> rVectors(nr*3);
	for (int ir = 0 ; ir < nr ; ++ir)
	{
		auto rVec = dvscf->get_real_space_grid().get_vector_direct(ir);
		for (int i = 0 ; i < 3 ; ++i)
			rVectors[ir*3+i] = rVec[i];
	}

	nM_ = ph->get_num_modes();
	nB_ = bandList.size();
	nBp_ = bandpList.size();

	Algorithms::LinearAlgebraInterface linalg;

	std::vector<double> allQVectors(nK_*nKp_*3);
	std::vector<std::pair<int,int>> q_to_k_and_kp_map(nK_*nKp_);
	for ( int ik = 0 ; ik < nK_; ++ik)
		for ( int ikp = 0 ; ikp < nKp_; ++ikp)
		{
			const int ikkp = ik*nKp_+ikp;
			allQVectors[ikkp*3+0] = kpList[ik*3+0]-kList[ikp*3+0];
			allQVectors[ikkp*3+1] = kpList[ik*3+1]-kList[ikp*3+1];
			allQVectors[ikkp*3+2] = kpList[ik*3+2]-kList[ikp*3+2];

			q_to_k_and_kp_map[ikkp].first = ik;
			q_to_k_and_kp_map[ikkp].second = ikp;
		}

	// map q vectors back to the 1. BZ
	for ( auto &qi : allQVectors )
		qi -= std::floor(qi+0.5);

	std::map<LatticeStructure::RegularBareGrid::GridPoint, std::vector<int>> gp_to_q_index;
	dvscf->query_q(allQVectors, gp_to_q_index);

	double percentageDone = 0.0;
	auto start_clock = std::chrono::system_clock::now();

	// for OpenMP, we need random access iterators ...
	std::vector<std::pair<LatticeStructure::RegularBareGrid::GridPoint, std::vector<int>>> gp_to_q_index_vector;
	gp_to_q_index_vector.reserve(gp_to_q_index.size());
	for ( auto const & coarseGP : gp_to_q_index )
		gp_to_q_index_vector.push_back(coarseGP);

#ifndef NDEBUG
	// cross check that all q-indices are there once.
	std::set<int> cross_check_set;
	for (int igp = 0 ; igp < gp_to_q_index_vector.size(); ++igp)
		for (int iq : gp_to_q_index_vector[igp].second )
			cross_check_set.insert(iq);
	assert(cross_check_set.size() == nK_*nKp_);
	assert(*cross_check_set.rbegin() == (nK_*nKp_-1));
#endif


	Algorithms::FFTInterface fft1, fft2;
	fft1.plan_fft(potentialFFTGrid, nB_,  1, false, nK_);
	fft2.plan_fft(potentialFFTGrid, nBp_, -1, false, nKp_);

	std::vector<int> localQDone(omp_get_max_threads(), 0);

	phononFrequencies_.resize(nK_*nKp_*nM_);
	data_.resize( nK_*nKp_*nM_*nB_*nBp_ );
	#pragma omp parallel
	{
		// thread private buffer data
		Auxillary::alignedvector::DV modes;
		Auxillary::alignedvector::ZV dynmat;
		Auxillary::alignedvector::CV wfcProdBuffRealSpace(nr), localGkkp, dvscfData, bufferWfct1,bufferWfct2;
		std::vector<Auxillary::alignedvector::CV> dvscfBuffers;

		#pragma omp for
		for (int icgp = 0 ; icgp < gp_to_q_index_vector.size(); ++icgp)
		{
		    int thread_id = omp_get_thread_num();

			auto const & coarseGP = gp_to_q_index_vector[icgp];

			// set things that go onto the coarse grid
			std::vector<double> const & gridQ = coarseGP.first.get_coords();

			ph->compute_at_q( gridQ, modes, dynmat );

			dvscf->compute_dvscf_q( gridQ, modes, dynmat, ph->get_masses(), dvscfData, dvscfBuffers);

			// now compute things that are handled without a grid
			for (int iq  : coarseGP.second )
			{
				int ik = q_to_k_and_kp_map[iq].first;
				int ikp = q_to_k_and_kp_map[iq].second;

				fft1.fft_sparse_data(
						fftMapsK[ik],
						wfcts->get_max_fft_dims(),
						wfcsk[ik],
						1,
						bufferWfct1);

				fft2.fft_sparse_data(
						fftMapsKp[ikp],
						wfcts->get_max_fft_dims(),
						wfcskp[ikp],
						-1,
						bufferWfct2);

				ph->compute_at_q( std::vector<double>(&allQVectors[iq*3], &allQVectors[iq*3]+3), modes, dynmat );
				std::copy(modes.begin(), modes.end(), &phononFrequencies_[(ik*nKp_+ikp)*nM_]);

				this->compute_gkkp_local(nr, nB_, nBp_, nM_, dvscfData, bufferWfct1, bufferWfct2, wfcProdBuffRealSpace, localGkkp);
				std::copy(localGkkp.begin(), localGkkp.end(), &data_[this->tensor_layout(ik,ikp,0,0,0)]);
			}

			// check if we crossed a cumulative number of njump q points
			const int njump = 1000;
			bool report = ((localQDone[thread_id] + coarseGP.second.size())%njump) < (localQDone[thread_id]%njump);

			localQDone[thread_id] += coarseGP.second.size();

			if (report and report_timings and (thread_id == 0))
			{
				percentageDone = 0.0;
				for (auto n : localQDone)
					percentageDone += static_cast<double>(n)/static_cast<double>(nK_*nKp_);
				auto thisTime_clock = std::chrono::system_clock::now();
				std::chrono::duration<double> elapsed_seconds = thisTime_clock-start_clock;
				std::chrono::duration<double> projectedRuntime_seconds = elapsed_seconds/percentageDone;
				std::chrono::duration<double> remainingRuntine = projectedRuntime_seconds-elapsed_seconds;

				std::cout << "\r\tPercentage done : " << std::floor(percentageDone*100 + 0.5)
						<< "; estimated time to finish: " << std::floor(remainingRuntine.count()+0.5) << " s"
						<<  "                 ";
				std::cout.flush();
			}
		}
	}//end #pragma omp parallel

	if ( report_timings )
	{
		std::cout << std::endl;
		for (int it = 0 ; it < localQDone.size() ; ++it)
			std::cout << "Thread: " << it << " finished " << localQDone[it] << " q-points" << std::endl;
	}
}

void
ElectronPhononCoupling::generate_gkkp_energy_units(
		std::vector<double> const & kList,
		std::vector<double> const & kpList,
		std::vector<int> bandList,
		std::vector<int> bandpList,
		std::shared_ptr<const Phonon> ph,
		std::shared_ptr<const DisplacementPotential> dvscf,
		std::shared_ptr<const ElectronicStructure::Wavefunctions> wfcts,
		std::vector<std::complex<float>> & gkkp) const
{
	assert(kList.size() == kpList.size());
	gkkp.clear();
	assert(ph and dvscf and wfcts);

	//Wave functions are on a regular grid and in G (reciprocal) space
	std::vector< std::vector<std::complex<float> > > wfcsk, wfcskp;
	std::vector< std::vector<int> > fftMapsK, fftMapsKp;
	wfcts->generate_wfcts_at_arbitray_kp( kList, bandList, wfcsk, fftMapsK);
	wfcts->generate_wfcts_at_arbitray_kp( kpList, bandpList, wfcskp, fftMapsKp );

	// conjugate the left wfct of the scalar product
	for ( auto & wfk : wfcsk )
		for ( auto & cg : wfk )
			cg = std::conj(cg);

	// some abbriviations
	std::vector<int> potentialFFTGrid = dvscf->get_real_space_grid().get_grid_dim();
	const int nk = kList.size()/3;
	const int nr = potentialFFTGrid[0]*potentialFFTGrid[1]*potentialFFTGrid[2];
	const int nB = bandList.size();
	const int nBp = bandpList.size();
	const int nM = dvscf->get_num_modes();

	std::vector<double> rVectors(nr*3);
	for (int ir = 0 ; ir < nr ; ++ir)
	{
		auto rVec = dvscf->get_real_space_grid().get_vector_direct(ir);
		for (int i = 0 ; i < 3 ; ++i)
			rVectors[ir*3+i] = rVec[i];
	}

	Algorithms::FFTInterface fft1, fft2;
	fft1.plan_fft(potentialFFTGrid, nB, 1, false, nk);
	fft2.plan_fft(potentialFFTGrid, nBp, -1, false, nk);

	std::vector<float> buffer(nr);
	Auxillary::alignedvector::DV modes;
	Auxillary::alignedvector::ZV dynmat;
	Auxillary::alignedvector::CV wfcProdBuffRealSpace(nr), localGkkp, dvscfData, bufferWfct1,bufferWfct2;
	std::vector<Auxillary::alignedvector::CV> dvscfBuffers;

	gkkp.reserve(nk*nB*nBp*nM);
	for ( int ik = 0 ; ik < nk ; ++ik )
	{
		// compute q vector in the 1. BZ
		std::vector<double> q{	kpList[ik*3+0]-kList[ik*3+0],
								kpList[ik*3+1]-kList[ik*3+1],
								kpList[ik*3+2]-kList[ik*3+2]};
		for ( auto &qi : q )
			qi -= std::floor(qi+0.5);

		// compute phonon modes and dynamical matrix
		ph->compute_at_q( q, modes, dynmat );
		dvscf->compute_dvscf_q(q, modes, dynmat, ph->get_masses(), dvscfData, dvscfBuffers);

		assert(modes.size() == nM);
		assert(dvscfData.size() == nr*nM);

		// generate wfcts on the potential grid
		fft1.fft_sparse_data(
				fftMapsK[ik],
				wfcts->get_max_fft_dims(),
				wfcsk[ik],
				1, // conjugate wfct
				bufferWfct1);

		fft2.fft_sparse_data(
				fftMapsKp[ik],
				wfcts->get_max_fft_dims(),
				wfcskp[ik],
				-1,
				bufferWfct2);

		this->compute_gkkp_local(nr, nB, nBp, nM, dvscfData, bufferWfct1, bufferWfct2, wfcProdBuffRealSpace, localGkkp);
		gkkp.insert(std::end(gkkp), std::begin(localGkkp), std::end(localGkkp));
	}
}

void
ElectronPhononCoupling::compute_gkkp_local(
		int nr,
		int nB,
		int nBp,
		int nM,
		Auxillary::alignedvector::CV const & dvscfData,
		Auxillary::alignedvector::CV const & wfctBufferk,
		Auxillary::alignedvector::CV const & wfctBufferkp,
		Auxillary::alignedvector::CV & wfcProdBuffRealSpace,
		Auxillary::alignedvector::CV & localGkkp) const
{
	assert( wfctBufferk.size() == nr*nB );
	assert( wfctBufferkp.size() == nr*nBp );
	assert( dvscfData.size() == nr*nM );
	assert( wfcProdBuffRealSpace.size() == nr);

	Algorithms::LinearAlgebraInterface linalg;
	localGkkp.assign(nB*nBp*nM, std::complex<float>(0.0f));

	for ( int ib = 0 ; ib < nB; ++ib )
		for ( int ibp = 0 ; ibp < nBp; ++ibp )
		{
			auto ptr_wf1 = wfctBufferk.data() + ib*nr;
			auto ptr_wf2 = wfctBufferkp.data() + ibp*nr;

			for (int ir = 0 ; ir < nr ; ++ir)
				wfcProdBuffRealSpace[ir] = ptr_wf1[ir]*ptr_wf2[ir]
											/static_cast<float>(nr);

			linalg.call_gemv('n', nM, nr,
					std::complex<float>(1.0f),
					dvscfData.data(), nM,
					wfcProdBuffRealSpace.data(), 1,
					std::complex<float>(0.0f),
					&localGkkp[this->local_tensor_layout(ib, ibp, 0, nB, nBp, nM)], 1);

//			for (int inu = 0 ; inu < nM; ++inu)
//				for (int ir = 0 ; ir < nr ; ++ir)
//					localGkkp[this->local_tensor_layout(ib, ibp, inu, nB, nBp, nM)] +=
//							dvscfData[inu*nM+ir]*wfcProdBuffRealSpace[ir];
		}
}

void
ElectronPhononCoupling::get_local_matrix_range(int ik, int ikp,
		std::vector< std::complex<float> >::iterator & rangeBegin,
		std::vector< std::complex<float> >::iterator & rangeEnd,
		std::vector< float >::iterator & phononFreqBegin,
		std::vector< float >::iterator & phononFreqEnd )
{
	rangeBegin = data_.begin()+this->tensor_layout(ik,ikp,0,0,0);
	rangeEnd = rangeBegin+nB_*nBp_*nM_;
	phononFreqBegin = phononFrequencies_.begin() + nM_*(ikp + nKp_*ik);
	phononFreqEnd = phononFrequencies_.begin() + nM_*(ikp + nKp_*ik) + nM_;
}

std::complex<float>
ElectronPhononCoupling::operator() (int ik, int ikp, int ib, int ibp, int imu) const
{
	int i = this->tensor_layout(ik, ikp, ib, ibp, imu);
	assert((data_.size()>i) && (i>=0));
	return data_[i];
}

int
ElectronPhononCoupling::tensor_layout(int ik, int ikp, int ib, int ibp, int imu) const
{
	assert( ik < nK_);
	assert( ikp < nKp_);
	assert( ib < nB_);
	assert( ibp < nBp_);
	assert( imu < nM_);
	return (((ik*nKp_+ikp)*nB_+ib)*nBp_+ibp)*nM_+imu;
}

int
ElectronPhononCoupling::local_tensor_layout( int ib, int ibp, int imu, int nB, int nBp, int nM) const
{
	assert( ib < nB);
	assert( ibp < nBp);
	assert( imu < nM);
	return (ib*nBp+ibp)*nM+imu;
}

void
ElectronPhononCoupling::write_gkkp_file(
		std::string const & filename,
		std::vector<double> const & kList,
		std::vector<double> const & kpList,
		std::vector<int> const & bandList,
		std::vector<int> const & bandpList) const
{
	std::ofstream file( filename.c_str(), std::ios::binary | std::ios::out );
	if ( not file.good() )
		throw std::runtime_error(std::string("Unable to open file ")+filename+" for writing the gkkp file.");

	// write the header
}


} /* namespace PhononStructure */
} /* namespace elephon */
