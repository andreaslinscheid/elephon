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

namespace elephon
{
namespace PhononStructure
{

void
ElectronPhononCoupling::generate_gkkp(
		std::vector<double> const & kList,
		std::vector<double> const & kpList,
		std::vector<int> const & bandList,
		std::vector<int> const & bandpList,
		Phonon const & ph,
		DisplacementPotential const & dvscf,
		ElectronicStructure::Wavefunctions const & wfcts)
{
	assert(kList.size()%3 == 0);
	assert(kpList.size()%3 == 0);
	nK_ = kList.size()/3;
	nKp_ = kpList.size()/3;

	std::vector<double> modes;
	std::vector<std::complex<double> > dynmat;
	std::vector<std::complex<float> > dvscfData;

	//Wave functions are on a regular grid and in G (reciprocal) space
	std::vector< std::vector<std::complex<float> > > wfcsk, wfcskp;
	std::vector< std::vector<int> > fftMapsK, fftMapsKp;
	wfcts.generate_wfcts_at_arbitray_kp( kList, bandList, wfcsk, fftMapsK);
	wfcts.generate_wfcts_at_arbitray_kp( kpList, bandpList, wfcskp, fftMapsKp );

	std::vector<int> potentialFFTGrid = dvscf.get_real_space_grid().get_grid_dim();

	std::vector<std::complex<float> > bufferWfct1,bufferWfct2;

	int nr = potentialFFTGrid[0]*potentialFFTGrid[1]*potentialFFTGrid[2];
	nM_ = ph.get_num_modes();
	nB_ = bandList.size();
	nBp_ = bandpList.size();

	data_.resize( nK_*nKp_*nM_*nB_*nBp_ );
	Algorithms::FFTInterface fft;
	for ( int ik = 0 ; ik < nK_; ++ik)
	{
		fft.fft_sparse_data(
				fftMapsK[ik],
				wfcts.get_max_fft_dims(),
				wfcsk[ik],
				nB_,
				-1,
				bufferWfct1,
				potentialFFTGrid,
				false,
				nK_*nKp_);

		for ( int ikp = 0 ; ikp < nKp_; ++ikp)
		{
			std::vector<double> q{	kList[ik*3+0]-kpList[ik*3+0],
									kList[ik*3+1]-kpList[ik*3+1],
									kList[ik*3+2]-kpList[ik*3+2]};
			ph.compute_at_q( q, modes, dynmat );
			dvscf.compute_dvscf_q( q, dynmat, ph.get_masses(), dvscfData);

			fft.fft_sparse_data(
					fftMapsKp[ik],
					wfcts.get_max_fft_dims(),
					wfcskp[ik],
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

					for (int imu = 0 ; imu < nM_ ; ++imu)
					{
						auto ptr_dvscf = dvscfData.data() + imu*nr;
						data_[this->tensor_layout(ik,ikp,ib,ibp,imu)] = std::complex<float>(0);
						for (int ir = 0 ; ir < nr ; ++ir)
							data_[this->tensor_layout(ik,ikp,ib,ibp,imu)] +=
									ptr_wf1[ir]*ptr_dvscf[ir]*ptr_wf2[ir];
					}
				}
		}
	}
}

void
ElectronPhononCoupling::get_local_matrix_range(int ik, int ikp,
		std::vector< std::complex<float> >::iterator & rangeBegin,
		std::vector< std::complex<float> >::iterator & rangeEnd )
{
	rangeBegin = data_.begin()+this->tensor_layout(ik,ikp,0,0,0);
	rangeEnd = data_.begin()+nB_*nBp_*nM_;
	assert( (data_.begin()+this->tensor_layout(ik,ikp+1,0,0,0)) == rangeEnd );
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

} /* namespace PhononStructure */
} /* namespace elephon */