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
	 * Compute all electron-phonon matrix elements squared between the vectors/band pairs k,i and k',i' for all modes.
	 *
	 * The units of the data is eV^2
	 *
	 * @param qVectorIndex	The regular grid index of the q vector
	 * @param reducibleIndicesKandKp
	 * @param bandList		List of bands i
	 * @param bandpList		List of bands i'
	 * @param ph			Pointer to the phonon module that allows to diagonalize the Fourier transformed matrix of force constants
	 * 						at every q = k'-k also yielding the dynamical matrix. Must be initialized.
	 * @param dvscf			Pointer to the displacement potential in real space that allows to compute the overlap
	 * 						created by the Fourier transformed displacement potential at q = k' - k between states |k,i> and |k',i'>
	 * @param wfcts			Pointer to the wave function module that allows to load and interpolate the wavefunction from a regular grid
	 * 						to a given non-grid point k for a band i.
	 * @param gkkpMod2		resized and set to the electron phonon matrix elements abs squared |<ki|dvscf(mu)|k'i'>|^2 with a layout
	 * 						where k is running slowest, then for each k the requested band matrix elements (ibnd,ibndp) and finally the modes.
	 */
	void generate_gkkp_mod_2_of_q(
			int qVectorIndex,
			std::vector<std::pair<int,int>> const & reducibleIndicesKandKp,
			std::vector<int> bandList,
			std::vector<int> bandpList,
			std::shared_ptr<const Phonon> ph,
			std::shared_ptr<const DisplacementPotential> dvscf,
			std::shared_ptr<const ElectronicStructure::Wavefunctions> wfcts,
			Auxillary::alignedvector::FV & gkkpMod2) const;

	void write_gkkp_file(
			std::string const & filename,
			std::vector<double> const & kList,
			std::vector<double> const & kpList,
			std::vector<int> const & bandList,
			std::vector<int> const & bandpList) const;

	/**
	 * For any give k and q points, this function defined the memory layout of the gkkp matrix produced by this class.
	 *
	 * @param ib	The index of the first band i in <ki|dvscf(mu)|k+qi'>
	 * @param ibp	The index of the second band i' in <ki|dvscf(mu)|k+qi'>
	 * @param imu	The index of the mode mu in <ki|dvscf(mu)|k+qi'>
	 * @param nB	The number of Bands i
	 * @param nBp	The number of Bands i'
	 * @param nM	The number of modes
	 * @return		The index of the element in the tensor <ki|dvscf(mu)|k+qi'> for given k and q.
	 */
	int local_tensor_layout(
			int ib,
			int ibp,
			int imu,
			int nB,
			int nBp,
			int nM) const;

private:

	void compute_gkkp_local(
			int nr,
			int nB,
			int nBp,
			int nM,
			Auxillary::alignedvector::CV const & dvscfData,
			Auxillary::alignedvector::CV const & wfctBufferk,
			Auxillary::alignedvector::CV const & wfctBufferkp,
			Auxillary::alignedvector::CV const & expGdot_r,
			Auxillary::alignedvector::CV & prodBuffer,
			Auxillary::alignedvector::CV & localGkkp) const;

	void compute_umklapp_phase(
			double gridPrec,
			std::vector<double> const &gVectors,
			LatticeStructure::RegularBareGrid const & rsGrid,
			std::vector<int> & buffMap,
			std::map<int,Auxillary::alignedvector::CV> & expGdot_r_buffer) const;
};

} /* namespace PhononStructure */
} /* namespace elephon */

#endif /* ELEPHON_PHONONSTRUCTURE_ELECTRONPHONONCOUPLING_H_ */
