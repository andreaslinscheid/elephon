/*	This file MatsubaraBaseFermion.h is part of elephon.
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
 *  Created on: May 16, 2018
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_ELIASHBERGEQUATIONS_MATSUBARABASEFERMION_H_
#define ELEPHON_ELIASHBERGEQUATIONS_MATSUBARABASEFERMION_H_

#include "EliashbergEquations/EliashbergModule.h"
#include "Auxillary/AlignedVector.h"

namespace elephon {
namespace EliashbergEquations {

/**
 * Base class for a Fermionic self-energy part.
 */
class MatsubaraBaseFermion
{
	typedef EliashbergModule::EliashbergDataType T;
public:

	/**
	 * Provide read/write access to the underlying container.
	 *
	 * The layout of the underlying data is the Matsubara index as the fast running index.
	 *
	 * @return ptr to the underlying data.
	 */
	T * access_data_ptr();

	/**
	 * Provide read access to the underlying container.
	 *
	 * The layout of the underlying data is the Matsubara index as the fast running index.
	 *
	 * @return constant ptr to the underlying data.
	 */
	T const * get_data_ptr() const;

	/**
	 * Pointer 1 past the end of the underlying data
	 *
	 * @return constant ptr to the underlying data.
	 */
	T const * get_end_data_ptr() const;

	/**
	 * Get the value at Matsubara frequency n and band b.
	 * @param[in] n		Matsubara frequency index
	 * @param[in] b		Band index.
	 * @return The value at this indices.
	 */
	T const & operator() (int n, int b) const;

	/**
	 * Set the value at Matsubara frequency n and band b.
	 * @param[in] n		Matsubara frequency index
	 * @param[in] b		Band index.
	 * @return A refernce to the value at this indices.
	 */
	T & operator() (int n, int b);

	/**
	 * Get the number of bands.
	 * @return the number of bands.
	 */
	int get_num_bands() const;

	/**
	 * Get the number of Matsubara frequencies for the current cutoff.
	 *
	 * @return the number of matsubara frequencies.
	 */
	int get_num_mats() const;

	/**
	 * Perform a spline interpolation to a new set of Matsubara frequencies, i.e. change the temperature.
	 *
	 * @param[in] oldMatsubaraFrequencies	Vector with the old Matsubara frequencies. Note that this object does
	 * 										not store the temperature explicitly.
	 * @param[in] newMatsubaraFrequencies	Vector with the new Matsubara frequencies
	 */
	void change_temperature(
			std::vector<double> const & oldMatsubaraFrequencies,
			std::vector<double> const & newMatsubaraFrequencies,
			std::shared_ptr<Auxillary::alignedvector::aligned_vector<double>> splineMatrixPtr = nullptr);
protected:

	/**
	 * Set the internal data.
	 * @param[in] numMatsubara	Number of Matsubara frequencies in the range of the energy cutoff
	 * @param[in] numBands		Number of bands
	 * @param[in] data			Vector that will be moved into the interal storage.
	 */
	void initialize(
			int numMatsubara,
			int numBands,
			Auxillary::alignedvector::aligned_vector<T> data);

private:
	int nMats_ = 0;

	int nDataPerMats_ = 0;

	Auxillary::alignedvector::aligned_vector<T> data_;
};

} /* namespace EliashbergEquations */
} /* namespace elephon */

#endif /* ELEPHON_ELIASHBERGEQUATIONS_MATSUBARABASEFERMION_H_ */
