/*	This file ASSymmetry.h is part of elephon.
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

#ifndef ELEPHON_ATOMICSITE_ASSYMMETRY_H_
#define ELEPHON_ATOMICSITE_ASSYMMETRY_H_

#include "AtomicSite/WignerDMatrix.h"
#include <vector>
#include <memory>

namespace elephon
{
namespace AtomicSite
{

/**
 *	Representation of the symmetry operations on the radial grid.
 */
class ASSymmetry
{
public:

	/**
	 * Type of a symmetry operation on a radial expansion.
	 */
	typedef std::vector<AtomicSite::WignerDMatrix> RadSym;

	/**
	 * Set the symmetry operations for up to within lmax based on the carthesian rotation matrices.
	 *
	 * @param[in] lmax							Maximal angular moment that will be considered. Thus we have a range of [0,lmax] rotation matrices.
	 * @param[in] carthesianSymmetryOperations	List of rotation matrices describing rotations in the passive sense (rotation of the coordinate system).
	 * 											Each block of 9 elements is interpreted as one C-layout matrix.
	 */
	void initialize(
			int lmax,
			std::vector<double> carthesianSymmetryOperations);

	/**
	 * Obtain a pointer to the set of Wigner D matrices for a symmetry operation.
	 *
	 * @param isymop	Index of the symmetry operation.
	 * @return			shared pointer to a constant vector of size lMax+1
	 * 						with rotation matrices from l=0 to l=lMax
	 */
	std::shared_ptr<const RadSym>
	get_wigner_rotation_matrices_symop(
			int isymop) const;

private:

	int lMax_ = 0;

	int numSymOps_ = 0;

	/// For each symmetry, the 3 Euler angles of the associated rotation
	/// The convention is z-y-z, as used for the Wigner D matrix
	/// See https://www.geometrictools.com/Documentation/EulerAngles.pdf
	/// and https://en.wikipedia.org/wiki/Wigner_D-matrix
	std::vector<double> eulerAngles_;

	/// For each symmetry operation, the representation in the space of
	/// spherical harmonics.
	std::vector<std::shared_ptr<RadSym>> rotMatricesPtr_;
};

} /* namespace AtomicSite */
} /* namespace elephon */

#endif /* ELEPHON_ATOMICSITE_ASSYMMETRY_H_ */
