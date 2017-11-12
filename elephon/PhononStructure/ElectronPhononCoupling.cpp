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
//	bool report_timings = false;

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

	std::vector<double> modes;
	std::vector<std::complex<double> > dynmat;
	std::vector<std::complex<float> > dvscfData;

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
		std::cout << "Generating wavefunctions at "<<kList.size()/3+kpList.size()/3
				<< " k points: "<<durations["wfcts_gen"].second.count()<<"s"<< std::endl;
	}

	std::vector<int> potentialFFTGrid = dvscf->get_real_space_grid().get_grid_dim();

	std::vector<std::complex<float> > bufferWfct1,bufferWfct2;

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

	std::vector<float> dotProdQr(nr), phasesRe(nr), phasesIm(nr), buffer(nr);
	std::vector<std::complex<float> > wfcProdBuffRealSpace(nr);
	std::vector<std::complex<float> > localGkkp;
	Algorithms::LinearAlgebraInterface linalg;

	phononFrequencies_.reserve(nK_*nKp_*nM_);
	data_.resize( nK_*nKp_*nM_*nB_*nBp_ );
	Algorithms::FFTInterface fft;
	auto start_clock = std::chrono::system_clock::now();
	for ( int ik = 0 ; ik < nK_; ++ik)
	{
		if ( report_timings )
			start_m_c("sparse_fft");
		fft.fft_sparse_data(
				fftMapsK[ik],
				wfcts->get_max_fft_dims(),
				wfcsk[ik],
				nB_,
				1,
				bufferWfct1,
				potentialFFTGrid,
				false,
				nK_*nKp_);

		if ( report_timings )
			stop_m_c("sparse_fft");

		for ( int ikp = 0 ; ikp < nKp_; ++ikp)
		{
			std::vector<double> q{	kList[ik*3+0]-kpList[ikp*3+0],
									kList[ik*3+1]-kpList[ikp*3+1],
									kList[ik*3+2]-kpList[ikp*3+2]};

			// map q vectors back to the 1. BZ
			for ( auto &qi : q )
				qi -= std::floor(qi+0.5);

			if ( report_timings )
				start_m_c("ph");
			ph->compute_at_q( q, modes, dynmat );
			if ( report_timings )
				stop_m_c("ph");
			if ( report_timings )
				start_m_c("ph_dvscf");
			dvscf->compute_dvscf_q( q, dynmat, ph->get_masses(), dvscfData);
			phononFrequencies_.insert(std::end(phononFrequencies_), modes.begin(), modes.end() );
			if ( report_timings )
				stop_m_c("ph_dvscf");

			if ( report_timings )
				start_m_c("sparse_fft");
			fft.fft_sparse_data(
					fftMapsKp[ikp],
					wfcts->get_max_fft_dims(),
					wfcskp[ikp],
					nBp_,
					-1,
					bufferWfct2,
					potentialFFTGrid,
					false,
					nK_*nKp_);
			if ( report_timings )
				stop_m_c("sparse_fft");

			// compute the phases from the Bloch phase factors of the two wave functions
			// NOTE: because this is the most performance critical part of the code, we
			//			jump throw some hoops to get it faster ...
			if ( report_timings )
				start_m_c("phases");
			this->compute_phases(nr, rVectors, q, buffer, phasesRe, phasesIm);

			if ( report_timings )
				stop_m_c("phases");

			if ( report_timings )
				start_m_c("gkkp_local");
			this->compute_gkkp_local(nr, nB_, nBp_, nM_, modes, dvscfData, bufferWfct1, bufferWfct2, phasesRe, phasesIm, localGkkp);
			if ( report_timings )
				stop_m_c("gkkp_local");

			std::copy(localGkkp.begin(), localGkkp.end(), &data_[this->tensor_layout(ik,ikp,0,0,0)]);
		}
		if ( report_timings )
		{
			auto thisTime_clock = std::chrono::system_clock::now();
			std::chrono::duration<double> elapsed_seconds = thisTime_clock-start_clock;
			std::chrono::duration<double> projectedRuntime_seconds = double(nK_) / double(ik+1) * elapsed_seconds;
			std::chrono::duration<double> remainingRuntine = projectedRuntime_seconds-elapsed_seconds;

			std::string timings;
			for ( auto &c : durations)
				timings += std::string(" ")+c.first + ": " + std::to_string(c.second.second.count()) + "; ";
			std::cout << "\r\tPercentage done : " << std::floor(100.0 *double(ik+1) / double(nK_) + 0.5)
					<< "; estimated time to finish: " << remainingRuntine.count() << "s\t"
					<< "timings: " << timings
					<<  "                 ";
			std::cout.flush();
		}
	}

	if ( report_timings )
		std::cout << std::endl;
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

	Algorithms::FFTInterface fft;

	std::vector<float> phasesRe(nr), phasesIm(nr), buffer(nr);
	std::vector<std::complex<float> > bufferWfct1,bufferWfct2;
	std::vector<double> modes;
	std::vector<std::complex<double> > dynmat;
	std::vector<std::complex<float> > dvscfData, localGkkp;
	gkkp.reserve(nk*nB*nBp*nM);
	for ( int ik = 0 ; ik < nk ; ++ik )
	{
		// compute q vector in the 1. BZ
		std::vector<double> q{	kList[ik*3+0]-kpList[ik*3+0],
								kList[ik*3+1]-kpList[ik*3+1],
								kList[ik*3+2]-kpList[ik*3+2]};
		for ( auto &qi : q )
			qi -= std::floor(qi+0.5);

		// compute the phases from the Bloch phase factors of the two wave functions
		this->compute_phases(nr, rVectors, q, buffer, phasesRe, phasesIm);

		// compute phonon modes and dynamical matrix
		ph->compute_at_q( q, modes, dynmat );
		dvscf->compute_dvscf_q( q, dynmat, ph->get_masses(), dvscfData);
		assert(modes.size() == nM);
		assert(dvscfData.size() == nr*nM);

		// generate wfcts on the potential grid
		fft.fft_sparse_data(
				fftMapsK[ik],
				wfcts->get_max_fft_dims(),
				wfcsk[ik],
				nB,
				1, // conjugate wfct
				bufferWfct1,
				potentialFFTGrid,
				/*dataLayoutRowMajor =*/false,
				nk*2);

		fft.fft_sparse_data(
				fftMapsKp[ik],
				wfcts->get_max_fft_dims(),
				wfcskp[ik],
				nBp,
				-1,
				bufferWfct2,
				potentialFFTGrid,
				/*dataLayoutRowMajor =*/false,
				nk*2);

		this->compute_gkkp_local(nr, nB, nBp, nM, modes, dvscfData, bufferWfct1, bufferWfct2, phasesRe, phasesIm, localGkkp);

		gkkp.insert(std::end(gkkp), std::begin(localGkkp), std::end(localGkkp));
	}
}

void
ElectronPhononCoupling::compute_phases(
		int nr,
		std::vector<double> const & rVectors,
		std::vector<double> const & q,
		std::vector<float> & buffer,
		std::vector<float> & phasesRe,
		std::vector<float> & phasesIm) const
{
	assert(buffer.size() == nr);
	assert(phasesRe.size() == nr);
	assert(phasesIm.size() == nr);
	assert(rVectors.size() == nr*3);
	assert(q.size() == 3);
	for ( int ir = 0 ; ir < nr ; ++ir)
		buffer[ir] = -2.0f*M_PI*(rVectors[ir*3+0]*q[0]+rVectors[ir*3+1]*q[1]+rVectors[ir*3+2]*q[2]);

	for ( int ir = 0 ; ir < nr ; ++ir)
	{
		phasesRe[ir] = std::cos(buffer[ir]);
		phasesIm[ir] = std::sin(buffer[ir]);
	}
}

void
ElectronPhononCoupling::compute_gkkp_local(
		int nr,
		int nB,
		int nBp,
		int nM,
		std::vector<double> const & modes,
		std::vector<std::complex<float>> const & dvscfData,
		std::vector<std::complex<float>> const & wfctBufferk,
		std::vector<std::complex<float>> const & wfctBufferkp,
		std::vector<float> const & phasesRe,
		std::vector<float> const & phasesIm,
		std::vector<std::complex<float>> & localGkkp) const
{
	assert( wfctBufferk.size() == nr*nB );
	assert( wfctBufferkp.size() == nr*nBp );
	assert( dvscfData.size() == nr*nM );
	assert( phasesIm.size() == nr);
	assert( phasesRe.size() == nr);

	Algorithms::LinearAlgebraInterface linalg;
	localGkkp.resize(nB*nBp*nM);
	std::vector<std::complex<float>> wfcProdBuffRealSpace(nr);

	for ( int ib = 0 ; ib < nB; ++ib )
		for ( int ibp = 0 ; ibp < nBp; ++ibp )
		{
			auto ptr_wf1 = wfctBufferk.data() + ib*nr;
			auto ptr_wf2 = wfctBufferkp.data() + ibp*nr;

			for (int ir = 0 ; ir < nr ; ++ir)
				wfcProdBuffRealSpace[ir] = ptr_wf1[ir]*ptr_wf2[ir]*std::complex<float>(phasesRe[ir], phasesIm[ir]);

			linalg.call_gemv('n', nM, nr,
					std::complex<float>(1.0f/nr),
					dvscfData.data(), nM,
					wfcProdBuffRealSpace.data(), 1,
					std::complex<float>(0.0f),
					&localGkkp[this->local_tensor_layout(ib, ibp, 0, nB, nBp, nM)], 1);

			// convert to displacement matrix element type, i.e. from units of eV/A to eV
			for (int inu = 0 ; inu < nM; ++inu)
			{
				double freq = modes[inu]; // modes is in units of Hz
				for (int ib = 0 ; ib < nB; ++ib)
					for (int ibp = 0 ; ibp < nBp; ++ibp)
						localGkkp[this->local_tensor_layout(ib, ibp, inu, nB, nBp, nM)] *=
									Auxillary::units::SQRT_HBAR_BY_2M_THZ_TO_ANGSTROEM/std::sqrt(freq);
			}
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
