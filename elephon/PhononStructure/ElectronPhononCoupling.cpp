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

	std::vector<double> modes;
	std::vector<std::complex<double> > dynmat;
	std::vector<std::complex<float> > dvscfData;

	//Wave functions are on a regular grid and in G (reciprocal) space
	std::vector< std::vector<std::complex<float> > > wfcsk, wfcskp;
	std::vector< std::vector<int> > fftMapsK, fftMapsKp;
	wfcts->generate_wfcts_at_arbitray_kp( kList, bandList, wfcsk, fftMapsK);
	wfcts->generate_wfcts_at_arbitray_kp( kpList, bandpList, wfcskp, fftMapsKp );

	// conjugate the left wfct of the scalar product
	for ( auto & wfk : wfcsk )
		for ( auto & cg : wfk )
			cg = std::conj(cg);

	std::vector<int> potentialFFTGrid = dvscf->get_real_space_grid().get_grid_dim();

	std::vector<std::complex<float> > bufferWfct1,bufferWfct2;

	int nr = potentialFFTGrid[0]*potentialFFTGrid[1]*potentialFFTGrid[2];
	nM_ = ph->get_num_modes();
	nB_ = bandList.size();
	nBp_ = bandpList.size();

	std::vector<std::complex<float> > wfcProdBuffRealSpace(nr);
	std::vector<std::complex<float> > localGkkp(nM_);
	Algorithms::LinearAlgebraInterface linalg;

	phononFrequencies_.reserve(nK_*nKp_*nM_);
	data_.resize( nK_*nKp_*nM_*nB_*nBp_ );
	Algorithms::FFTInterface fft;
	auto start_clock = std::chrono::system_clock::now();
	for ( int ik = 0 ; ik < nK_; ++ik)
	{
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

		for ( int ikp = 0 ; ikp < nKp_; ++ikp)
		{
			std::vector<double> q{	kList[ik*3+0]-kpList[ik*3+0],
									kList[ik*3+1]-kpList[ik*3+1],
									kList[ik*3+2]-kpList[ik*3+2]};

			// map q vectors back to the 1. BZ
			for ( auto &qi : q )
				qi -= std::floor(qi+0.5);

			ph->compute_at_q( q, modes, dynmat );
			dvscf->compute_dvscf_q( q, dynmat, ph->get_masses(), dvscfData);
			phononFrequencies_.insert(std::end(phononFrequencies_), modes.begin(), modes.end() );

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

			assert( bufferWfct1.size() == nr*nB_ );
			assert( bufferWfct2.size() == nr*nBp_ );
			assert( dvscfData.size() == nr*nM_ );

			for ( int ib = 0 ; ib < bandList.size(); ++ib )
				for ( int ibp = 0 ; ibp < bandpList.size(); ++ibp )
				{
					auto ptr_wf1 = bufferWfct1.data() + ib*nr;
					auto ptr_wf2 = bufferWfct2.data() + ibp*nr;

					for (int ir = 0 ; ir < nr ; ++ir)
						wfcProdBuffRealSpace[ir] = ptr_wf1[ir]*ptr_wf2[ir];

					linalg.call_gemv('n', nM_, nr,
							std::complex<float>(1.0f/static_cast<float>(nr)),
							dvscfData.data(), nM_,
							wfcProdBuffRealSpace.data(), 1,
							std::complex<float>(0.0f), localGkkp.data(), 1);
					for (int imu = 0 ; imu < nM_ ; ++imu)
						data_[this->tensor_layout(ik,ikp,ib,ibp,imu)] = localGkkp[imu];
				}
		}
		auto thisTime_clock = std::chrono::system_clock::now();
	    std::chrono::duration<double> elapsed_seconds = thisTime_clock-start_clock;
	    std::chrono::duration<double> projectedRuntime_seconds = double(nK_) / double(ik+1) * elapsed_seconds;
	    std::chrono::duration<double> remainingRuntine = projectedRuntime_seconds-elapsed_seconds;

		std::cout << "\r\tPercentage done : " << std::floor(100.0 *double(ik+1) / double(nK_) + 0.5)
				<< "; estimated time to finish: " << remainingRuntine.count() << "s"
				<<  "                 ";
		std::cout.flush();
	}
	std::cout << std::endl;
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
