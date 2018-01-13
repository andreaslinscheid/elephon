/*	This file SphericalHarmonicExpansion.h is part of elephon.
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
 *  Created on: Jan 4, 2018
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_ATOMICSITE_SPHERICALHARMONICEXPANSION_H_
#define ELEPHON_ATOMICSITE_SPHERICALHARMONICEXPANSION_H_

#include "Auxillary/AlignedVector.h"
#include "LatticeStructure/Symmetry.h"
#include "AtomicSite/RadialGrid.h"
#include "AtomicSite/WignerDMatrix.h"

namespace elephon
{
namespace AtomicSite
{

/**
 * A container for a spherical harmonic expansion.
 */
class SphericalHarmonicExpansion
{
public:

	/**
	 * Set the internal data.
	 *
	 * @param lmax	maximal angular momentum that occurs l=(0,...,lmax)
	 * @param rMax	dimension of the data per l and m
	 * @param data	list of data. Must be of exactly the size determined by lmax and rMax.
	 * 				The layout of the data must be r fastest and l slowest, such that
	 * 				for each constant l we have a block of m=-l to m=l which consists of a block of rMax elements,
	 * 				Since the sum of natural numbers is l(l+1)/2, we can compute the location of an element
	 * 				if r ranges from 0 to rMax, m from -l to l and l from 0 to lmax as
	 * 				r+rMax*( sum _{l' < l} (2l'+1) + m+l ) = r+rMax*(l*l + m+l)
	 */
	void initialize(
			int lmax,
			int rMax,
			Auxillary::alignedvector::ZV data,
			RadialGrid rgrid);

	/**
	 * Access an element of the expansion
	 *
	 * @param r		remainder index per l and m quantum number
	 * @param l		angular quantum number l of the spherical harmonic
	 * @param m		magnetic quantum number m of the spherical harmonic, range is m=[-l,l]
	 * @return		a reference to the element specified by the indices.
	 */
	std::complex<double> & operator() (int r, int m, int l);

	/**
	 * Obtain an element of the expansion.
	 *
	 * @param r		remainder index per l and m quantum number
	 * @param l		angular quantum number l of the spherical harmonic
	 * @param m		magnetic quantum number m of the spherical harmonic, range is m=[-l,l]
	 * @return		a copy of the element specified by the indices.
	 */
	std::complex<double> operator() (int r, int m, int l) const;

	/**
	 *	Defines the memory layout of angular momentum channels.
	 *
	 * @param l		The main angular momentum quantum number.
	 * @param m		The magnetic quantum number, m is in the range [-l,l]
	 * @return		The position of the element (l,m)
	 */
	int angular_momentum_layout(int l, int m) const;

	/**
	 * Perform a rotation of the data in this object according Wigner rotation matrices.
	 *
	 * @param wignerD	A vector with the rotation matrices for l=0,1,...,Lmax in this order.
	 * 					The maximal Lmax must be larger or equal than the internal range of the
	 * 					data in this object.
	 */
	void apply_wigner_D_rotation(std::vector<WignerDMatrix> const & wignerD);

	void transform(LatticeStructure::Symmetry::SymmetryOperation const & sop);

	void interpolate(
			std::vector<double> const & coordinates,
			Auxillary::alignedvector::ZV & interpolated_data) const;

private:

	int lmax_ = 0;

	int rMax_ = 0;

	Auxillary::alignedvector::ZV data_;

	RadialGrid rgrid_;
};

} /* namespace AtomicSite */
} /* namespace elephon */

#endif /* ELEPHON_ATOMICSITE_SPHERICALHARMONICEXPANSION_H_ */
