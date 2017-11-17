/*	This file ElectronPhononCouling.h is part of elephon.
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

#ifndef ELEPHON_PHONONSTRUCTURE_ELECTRONPHONONCOUPLING_H_
#define ELEPHON_PHONONSTRUCTURE_ELECTRONPHONONCOUPLING_H_

#include "PhononStructure/Phonon.h"
#include "PhononStructure/DisplacementPotential.h"
#include "ElectronicStructure/Wavefunctions.h"
#include "Auxillary/AlignedVector.h"
#include <vector>
#include <complex>
#include <memory>

namespace elephon
{
namespace PhononStructure
{

class ElectronPhononCoupling
{
public:

	/**
	 * Compute all electron-phonon matrix elements between all combination of vectors/band pairs k,i and k',i' for all modes.
	 *
	 * After this method was called, the internal storage is set to the complex gkkp matrix elements. For efficiency these are single
	 * precision complex.
	 * The units of the internal object set are eV/A.
	 *
	 * @param kList			List of k1x,k1y,k1z , k2x ... coordinates of the k vectors in units of the reciprocal lattice
	 * @param kpList		List of k1x,k1y,k1z , k2x ... coordinates of the k' vectors in units of the reciprocal lattice
	 * @param bandList		List of bands i
	 * @param bandpList		List of bands i'
	 * @param ph			Pointer to the phonon module that allows to diagonalize the Fourier transformed matrix of force constants
	 * 						at every q = k'-k also yielding the dynamical matrix. Must be initialized.
	 * @param dvscf			Pointer to the displacement potential in real space that allows to compute the overlap
	 * 						created by the Fourier transformed displacement potential at q = k' - k between states |k,i> and |k',i'>
	 * @param wfcts			Pointer to the wave function module that allows to load and interpolate the wavefunction from a regular grid
	 * 						to a given non-grid point k for a band i.
	 */
	void generate_gkkp_and_phonon(
			std::vector<double> kList,
			std::vector<double> kpList,
			std::vector<int> bandList,
			std::vector<int> bandpList,
			std::shared_ptr<const Phonon> ph,
			std::shared_ptr<const DisplacementPotential> dvscf,
			std::shared_ptr<const ElectronicStructure::Wavefunctions> wfcts );

	/**
	 * Compute selected electron-phonon matrix elements between the vectors and bands specified for all modes.
	 *
	 * This variant only computes <k1|dvscf|k1'>, <k2|dvscf|k2'>, ... so not all combinations of k and k' vectors.
	 *
	 * @param kList			List of k1x,k1y,k1z , k2x ... coordinates of the k vectors in units of the reciprocal lattice
	 * @param kpList		List of k1x,k1y,k1z , k2x ... coordinates of the k' vectors in units of the reciprocal lattice.
	 * 						Since exactly <ki|dvscf(mu)|k'i'> is computed, kpList must have the same length as \p kList.
	 * @param bandList		List of bands i
	 * @param bandpList		List of bands i'
	 * @param ph			Pointer to the phonon module that allows to diagonalize the Fourier transformed matrix of force constants
	 * 						at every q = k'-k also yielding the dynamical matrix. Must be initialized.
	 * @param dvscf			Pointer to the displacement potential in real space that allows to compute the overlap
	 * 						created by the Fourier transformed displacement potential at q = k' - k between states |k,i> and |k',i'>
	 * @param wfcts			Pointer to the wave function module that allows to load and interpolate the wavefunction from a regular grid
	 * 						to a given non-grid point k for a band i.
	 * @param gkkp			resized and set to the electron phonon matrix elements <ki|dvscf(mu)|k'i'> for each with a layout
	 * 						where k is running slowest, then for each k the request band matrix elements and finally the modes.
	 */
	void generate_gkkp_energy_units(
			std::vector<double> const & kList,
			std::vector<double> const & kpList,
			std::vector<int> bandList,
			std::vector<int> bandpList,
			std::shared_ptr<const Phonon> ph,
			std::shared_ptr<const DisplacementPotential> dvscf,
			std::shared_ptr<const ElectronicStructure::Wavefunctions> wfcts,
			std::vector<std::complex<float>> & gkkp) const;

	void write_gkkp_file(
			std::string const & filename,
			std::vector<double> const & kList,
			std::vector<double> const & kpList,
			std::vector<int> const & bandList,
			std::vector<int> const & bandpList) const;

	std::complex<float> operator() (int ik, int ikp, int ib, int ibp, int imu) const;

	void get_local_matrix_range(int ik, int ikp,
			std::vector< std::complex<float> >::iterator & rangeBegin,
			std::vector< std::complex<float> >::iterator & rangeEnd,
			std::vector< float >::iterator & phononFreqBegin,
			std::vector< float >::iterator & phononFreqEnd);

private:

	int nK_ = 0;

	int nB_ = 0;

	int nKp_ = 0;

	int nBp_ = 0;

	int nM_ = 0;

	std::vector<std::complex<float>> data_;

	std::vector<float> phononFrequencies_;

	int tensor_layout(int ik, int ikp, int ib, int ibp, int imu) const;

	int local_tensor_layout( int ib, int ibp, int imu, int nB, int nBp, int nM) const;

	void compute_gkkp_local(
			int nr,
			int nB,
			int nBp,
			int nM,
			Auxillary::alignedvector::CV const & dvscfData,
			Auxillary::alignedvector::CV const & wfctBufferk,
			Auxillary::alignedvector::CV const & wfctBufferkp,
			Auxillary::alignedvector::CV & wfcProdBuffRealSpace,
			Auxillary::alignedvector::CV & localGkkp) const;
};

} /* namespace PhononStructure */
} /* namespace elephon */

#endif /* ELEPHON_PHONONSTRUCTURE_ELECTRONPHONONCOUPLING_H_ */
