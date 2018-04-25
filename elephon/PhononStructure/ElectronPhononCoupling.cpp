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
ElectronPhononCoupling::generate_gkkp_mod_2_of_q(
		int qVectorIndex,
		std::vector<std::pair<int,int>> const & reducibleIndicesKandKp,
		std::vector<int> bandList,
		std::vector<int> bandpList,
		std::shared_ptr<const Phonon> ph,
		std::shared_ptr<const DisplacementPotential> dvscf,
		std::shared_ptr<const ElectronicStructure::Wavefunctions> wfcts,
		Auxillary::alignedvector::FV & gkkpMod2) const
{
	auto kGrid = wfcts->get_k_grid().view_bare_grid();
	const int nK = reducibleIndicesKandKp.size();
	const int nB = bandList.size();
	const int nBp = bandpList.size();
	const int nM = dvscf->get_num_modes();

#ifndef NDEBUG
	auto xyzQ = kGrid.get_reducible_to_xyz(qVectorIndex);
#endif
	std::vector<double> gridQ = kGrid.get_vector_direct(qVectorIndex);

	// create two vectors from the input k index pairs
	std::vector<int> reducibleIndicesK(nK);
	std::vector<int> reducibleIndicesKp(nK);
	for (int ikp = 0 ; ikp < nK; ++ikp)
	{
		reducibleIndicesK[ikp]=reducibleIndicesKandKp[ikp].first;
		reducibleIndicesKp[ikp]=reducibleIndicesKandKp[ikp].second;
	}

	std::vector<double> reducibleIndicesG(reducibleIndicesK.size()*3);
	std::vector<int> xyz(3);
	std::vector<double> v(3);
	// q = k' - k + G0
	for (int ikp = 0 ; ikp < reducibleIndicesK.size(); ++ikp)
	{
		// compute k'-k
		kGrid.get_vector_direct(reducibleIndicesKp[ikp], xyz, v);
		for (int i = 0 ; i < 3; ++i)
			reducibleIndicesG[ikp*3+i] = v[i];
		kGrid.get_vector_direct(reducibleIndicesK[ikp], xyz, v);
		for (int i = 0 ; i < 3; ++i)
			reducibleIndicesG[ikp*3+i] -= v[i];
#ifndef NDEBUG
		// cross check that k+q = k' modulo a grid vector.
		kGrid.get_reducible_to_xyz(reducibleIndicesK[ikp], xyz);
		for (int i = 0 ; i < 3; ++i)
			xyz[i] += xyzQ[i];
		const int ikcheck = kGrid.get_xyz_to_reducible_periodic(xyz);
		assert(reducibleIndicesKp[ikp] == ikcheck);
#endif
	}
	// now reducibleIndicesG holds q + G0. Get G0
	for( auto &gi : reducibleIndicesG )
		gi = std::floor(gi+0.5);

	//Wave functions are on a regular grid and in G (reciprocal) space
	//TODO: optimization opportunity: consider merging the loads for equivalent k points
	std::vector< std::vector<std::complex<float> > > wfcskp, wfcsk;
	std::vector<int> npwPerK, npwPerKp;
	wfcts->generate_reducible_grid_wfcts(bandList, reducibleIndicesK, wfcsk, npwPerK );
	wfcts->generate_reducible_grid_wfcts(bandpList, reducibleIndicesKp, wfcskp, npwPerKp);

	// conjugate the left wfct of the scalar product
	for ( auto & wfk : wfcsk )
		for ( auto & cg : wfk )
			cg = std::conj(cg);

	std::vector<double> kVectors = kGrid.get_vectors_direct(reducibleIndicesK);
	std::vector<double> kpVectors = kGrid.get_vectors_direct(reducibleIndicesKp);

	std::vector< std::vector<int> > fftMapsKp, fftMapsK;
	wfcts->compute_Fourier_maps(kpVectors, fftMapsKp);
	wfcts->compute_Fourier_maps(kVectors, fftMapsK);

	auto const & realSpaceGrid = dvscf->get_real_space_grid();
	std::vector<int> potentialFFTGrid = realSpaceGrid.get_grid_dim();

	int nr = potentialFFTGrid[0]*potentialFFTGrid[1]*potentialFFTGrid[2];

	Algorithms::FFTInterface fft1, fft2;
	fft1.plan_fft(potentialFFTGrid, nB,  1, false, nK);
	fft2.plan_fft(potentialFFTGrid, nBp, -1, false, nK);

	Auxillary::Multi_array<double,2> modes;
	Auxillary::Multi_array<std::complex<double>,3> dynmat, conjDynmat;
	Auxillary::alignedvector::CV localGkkp, dvscfData, bufferWfct1,bufferWfct2, phaseBuffer(nr);
	std::vector<Auxillary::alignedvector::CV> dvscfBuffers;
	Auxillary::alignedvector::ZV localGkkpMod2(nM*nM), localGkkpMod2Buffer(nM*nM);

	ph->compute_at_q( gridQ, modes, dynmat );
	dvscf->compute_dvscf_q( gridQ, modes, dynmat, ph->get_masses(), dvscfData, dvscfBuffers);

	// Pre-compute all the Umklapp phases that occur
	std::map<int, Auxillary::alignedvector::CV> expG0dot_r_buffer;
	std::vector<int> bufferMap;
	this->compute_umklapp_phase(
			kGrid.get_grid_prec(),
			reducibleIndicesG,
			realSpaceGrid,
			bufferMap,
			expG0dot_r_buffer);

	gkkpMod2.assign(nK*nB*nBp*nM, 0.0f);
	for ( int ik = 0 ; ik < nK; ++ik)
	{
		fft1.fft_sparse_data(
				fftMapsK[ik],
				wfcts->get_max_fft_dims(),
				wfcsk[ik],
				1,
				bufferWfct1);

		fft2.fft_sparse_data(
				fftMapsKp[ik],
				wfcts->get_max_fft_dims(),
				wfcskp[ik],
				-1,
				bufferWfct2);

		this->compute_gkkp_local(
				nr,
				nB,
				nBp,
				nM,
				dvscfData,
				bufferWfct1,
				bufferWfct2,
				expG0dot_r_buffer[bufferMap[ik]],
				phaseBuffer,
				localGkkp);

		for (int ib = 0 ; ib < nB ; ++ib)
			for (int ibp = 0 ; ibp < nBp ; ++ibp)
				for (int imu = 0 ; imu < nM ; ++imu)
				{
					gkkpMod2[((ik*nB+ib)*nBp+ibp)*nM+imu] =
							std::real(std::conj(localGkkp[(ib*nBp+ibp)*nM+imu])*localGkkp[(ib*nBp+ibp)*nM+imu]);
					assert(gkkpMod2[((ik*nB+ib)*nBp+ibp)*nM+imu] == gkkpMod2[((ik*nB+ib)*nBp+ibp)*nM+imu]);
				}

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
		Auxillary::alignedvector::CV const & expGdot_r,
		Auxillary::alignedvector::CV & prodBuffer,
		Auxillary::alignedvector::CV & localGkkp) const
{
	assert( wfctBufferk.size() == nr*nB );
	assert( wfctBufferkp.size() == nr*nBp );
	assert( dvscfData.size() == nr*nM );
	assert( expGdot_r.size() == nr);
	assert( prodBuffer.size() == nr);

	Algorithms::LinearAlgebraInterface linalg;
	localGkkp.assign(nB*nBp*nM, std::complex<float>(0.0f));

	for ( int ib = 0 ; ib < nB; ++ib )
		for ( int ibp = 0 ; ibp < nBp; ++ibp )
		{
			auto ptr_wf1 = wfctBufferk.data() + ib*nr;
			auto ptr_wf2 = wfctBufferkp.data() + ibp*nr;

			std::copy(expGdot_r.begin(), expGdot_r.end(), prodBuffer.begin());

			for (int ir = 0 ; ir < nr ; ++ir)
				prodBuffer[ir] *= ptr_wf1[ir]*ptr_wf2[ir]
											/static_cast<float>(nr);

			linalg.call_gemv('n', nM, nr,
					std::complex<float>(1.0f),
					dvscfData.data(), nr,
					prodBuffer.data(), 1,
					std::complex<float>(0.0f),
					&localGkkp[this->local_tensor_layout(ib, ibp, 0, nB, nBp, nM)], 1);
		}
}

int
ElectronPhononCoupling::local_tensor_layout(
		int ib,
		int ibp,
		int imu,
		int nB,
		int nBp,
		int nM) const
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

void
ElectronPhononCoupling::compute_umklapp_phase(
		double gridPrec,
		std::vector<double> const &gVectors,
		LatticeStructure::RegularBareGrid const & rsGrid,
		std::vector<int> & buffMap,
		std::map<int,Auxillary::alignedvector::CV> & expGdot_r_buffer) const
{
	const int nK =  gVectors.size()/3;
	buffMap.resize(nK);

	assert(*std::max_element(gVectors.begin(), gVectors.end()) < 1+gridPrec);
	assert(*std::min_element(gVectors.begin(), gVectors.end()) > -1-gridPrec);
	auto map_int = [&] (double val) {
		int ival = std::floor(val+0.5);
		return ival < 0 ? ival+3 : ival;};

	// the phase is trivially 1 if G is zero
	expGdot_r_buffer[0] = std::move(Auxillary::alignedvector::CV(rsGrid.get_num_points(),std::complex<float>(1.0f)));

	for (int ik = 0 ; ik < nK ; ++ik)
	{
		const int iGx =  map_int(gVectors[ik*3+0]);
		const int iGy =  map_int(gVectors[ik*3+1]);
		const int iGz =  map_int(gVectors[ik*3+2]);
		const int buffindex = iGx+2*(iGy+2*iGz);
		buffMap[ik] = buffindex;
		if ( expGdot_r_buffer.find(buffindex) == expGdot_r_buffer.end())
		{
			// the code below attempts to reduce expensive complex phase calculations in this performance critical part of the code
			const int nrx = rsGrid.get_grid_dim()[0];
			const int nry = rsGrid.get_grid_dim()[1];
			const int nrz = rsGrid.get_grid_dim()[2];
			auto phaseIncX = std::complex<double>( std::cos(2.0*M_PI*gVectors[ik*3+0]/static_cast<double>(nrx)),
												   std::sin(2.0*M_PI*gVectors[ik*3+0]/static_cast<double>(nrx)) );
			auto phaseIncY = std::complex<double>( std::cos(2.0*M_PI*gVectors[ik*3+1]/static_cast<double>(nry)),
												   std::sin(2.0*M_PI*gVectors[ik*3+1]/static_cast<double>(nry)) );
			auto phaseIncZ = std::complex<double>( std::cos(2.0*M_PI*gVectors[ik*3+2]/static_cast<double>(nrz)),
												   std::sin(2.0*M_PI*gVectors[ik*3+2]/static_cast<double>(nrz)) );
			Auxillary::alignedvector::CV phases(rsGrid.get_num_points());

			std::complex<double> totalPhaseGP(1.0);
			for (int iz = 0 ; iz < nrz; ++iz)
			{
				for (int iy = 0 ; iy < nry; ++iy)
				{
					for (int ix = 0 ; ix < nrx; ++ix)
					{
						const int ir = ix+nrx*(iy+nry*iz);
						phases[ir] = static_cast<std::complex<float>>(totalPhaseGP);
#ifndef NDEBUG
						// We use explicitly that the grid is x major. This assertion should fail if anybody ever was to change that
						auto rVec = rsGrid.get_vector_direct(ir);
						auto pha = std::exp(std::complex<double>(0.0, 2.0*M_PI*(gVectors[ik*3+0]*rVec[0]+
																				gVectors[ik*3+1]*rVec[1]+
																				gVectors[ik*3+2]*rVec[2] )));
						assert(std::abs(totalPhaseGP-pha)<1e-8);
#endif
						totalPhaseGP *= phaseIncX;
					}
					totalPhaseGP *= phaseIncY;
				}
				totalPhaseGP *= phaseIncZ;
			}
			expGdot_r_buffer[buffindex] = std::move(phases);
		}
	}
}

} /* namespace PhononStructure */
} /* namespace elephon */
