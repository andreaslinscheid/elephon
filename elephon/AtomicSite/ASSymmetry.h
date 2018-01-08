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

#include "LatticeStructure/LatticeModule.h"
#include "LatticeStructure/Atom.h"
#include "AtomicSite/WignerDMatrix.h"
#include <vector>
#include <memory>

namespace elephon
{
namespace AtomicSite
{

/**
 *
 */
class ASSymmetry
{
public:

	void initialize(
			int lmax,
			std::vector<int> symmetryOperations,
			std::vector<double> fractionalTranslations,
			LatticeStructure::LatticeModule lattice,
			std::vector<LatticeStructure::Atom> const & atoms);

	int get_num_sym_ops() const;

//	void transform(
//			int isym,
//			SphericalHarmonicExpansion & data) const;

	void get_euler_angles(
			int isym,
			double & alpha,
			double & beta,
			double & gamma) const;

	/**
	 * Obtain a pointer to the set of Wigner D matrices for a symmetry operation.
	 *
	 * @param isymop	Index of the symmetry operation.
	 * @return			shared pointer to a constant vector of size lMax+1
	 * 						with rotation matrices from l=0 to l=lMax
	 */
	std::shared_ptr<const std::vector<AtomicSite::WignerDMatrix>>
	get_wigner_rotation_matrices_symop(
			int isymop) const;

private:

	int lMax_ = 0;

	int numSymOps_ = 0;

	/// For each symmetry, the position in the list of
	/// atomic sites a given atom is moved to under application of this symmetry.
	/// Example: siteMapSymmetry_[1][0] == 1 means under symmetry op no. 1
	///          atom 0 is moved to the site of atom 1
	std::vector<std::vector<int>> siteMapSymmetry_;

	/// For each symmetry, the 3 Euler angles of the associated rotation
	/// The convention is z-y-z, as used for the Wigner D matrix
	/// See https://www.geometrictools.com/Documentation/EulerAngles.pdf
	/// and https://en.wikipedia.org/wiki/Wigner_D-matrix
	std::vector<double> eulerAngles_;

	std::vector<std::shared_ptr<std::vector<AtomicSite::WignerDMatrix>>> rotMatricesPtr_;
};

} /* namespace AtomicSite */
} /* namespace elephon */

#endif /* ELEPHON_ATOMICSITE_ASSYMMETRY_H_ */
