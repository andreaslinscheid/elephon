/*	This file WignerDMatrix.h is part of elephon.
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
 *  Created on: Jan 2, 2018
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_ATOMICSITE_WIGNERDMATRIX_H_
#define ELEPHON_ATOMICSITE_WIGNERDMATRIX_H_

#include "Auxillary/AlignedVector.h"
#include <vector>
#include <complex>

namespace elephon
{
namespace AtomicSite
{

/**
 * A rotation matrix for spherical harmonics.
 */
class WignerDMatrix
{
public:
	/**
	 * Fill the internal storage with the (2*l+1)x(2l+1) matrix elements of the Wigner D rotation matrix
	 *
	 * The convention is such that according to Wikipedia [https://en.wikipedia.org/wiki/Wigner_D-matrix]
	 * we are computing the complex conjugate. The convention means, we have the following transformation
	 * property for spherical harmonics:
	 *   Ylm(r') = sum m'=-l,...,+l D^l(alpha, beta, gamma)_m,m' Ylm'(r)
	 * l is the angular quantum number
	 *
	 * @param angularQuantumNumber	dimensionality of the representation, l=0 scalar, l=1 vector ...
	 * @param alpha					Euler angle of the second rotation about the z axis
	 * @param beta					Euler angle of the rotation about the y axis
	 * @param gamma					Euler angle of the first rotation about the z axis
	 */
	void initialize(
			int angularQuantumNumber,
			double alpha,
			double beta,
			double gamma);

	/**
	 * Obtain a single element of the Wigner matrix.
	 *
	 * @param m		Magnetic quantum number 1, range is [-l,l]
	 * @param mp	Magnetic quantum number 2, range is [-l,l]
	 * @return		A copy of the element D^(l)_{m,mp}(alpha, beta, gamma)
	 */
	std::complex<double> operator() (int m, int mp) const;

	/**
	 * View the internal (2l+1)x(2l+1) data.
	 *
	 * @return constant reference the internal (2l+1)x(2l+1) data
	 */
	Auxillary::alignedvector::ZV const &
	view_as_matrix() const;

private:

	/// The angular quantum number
	int l_ = 0;

	/// The Wigner D matrix has dimension 2l+1
	/// i.e. angularQuantumNumber*2 + 1
	int dim_ = 0;

	/// Storage of the full Wigner D matrix
	Auxillary::alignedvector::ZV matrix_;

	void wigner_small_d(
			double beta,
			std::vector<double> & D) const;

	void wigner_small_d_first_quad(
			double beta,
			std::vector<double> & D) const;

	int angular_layout(int m, int mp) const;

	double power_minus_1 (int m) const;
};

} /* namespace AtomicSite */
} /* namespace elephon */

#endif /* ELEPHON_ATOMICSITE_WIGNERDMATRIX_H_ */
