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
#include "AtomicSite/RadialGrid.h"

namespace elephon
{
// forward declares
namespace symmetry { class SymmetryOperation; };

namespace AtomicSite
{

/**
 * A container for a spherical harmonic expansion.
 *
 * It provides the basic interface to handle numerical data in this representation
 * as well as methods to interpolate back onto other points via SphericalHarmonicExpansion::interpolate
 * and fitting from a different data representation onto the spherical harmonic expansion via the method
 * SphericalHarmonicExpansion::fit_to_data.
 */
class SphericalHarmonicExpansion
{
public:

	/**
	 * Set the internal data.
	 *
	 * @param[in] lmax	maximal angular momentum that occurs l=(0,...,lmax)
	 * @param[in] data	list of data. Must be of exactly the size determined by lmax and rMax.
	 * 					The layout of the data must be r fastest and l slowest, such that
	 * 					for each constant l we have a block of m=-l to m=l which consists of a block of rMax elements,
	 * 					Since the sum of natural numbers is l(l+1)/2, we can compute the location of an element
	 * 					if r ranges from 0 to rMax, m from -l to l and l from 0 to lmax as
	 * 					r+rMax*( sum _{l' < l} (2l'+1) + m+l ) = r+rMax*(l*l + m+l)
	 * @param[in] rgrid	The radial grid.
	 */
	void initialize(
			int lmax,
			Auxillary::alignedvector::ZV data,
			RadialGrid rgrid);

	/**
	 * Access an element of the expansion
	 *
	 * @param[in] r		remainder index per l and m quantum number
	 * @param[in] l		angular quantum number l of the spherical harmonic
	 * @param[in] m		magnetic quantum number m of the spherical harmonic, range is m=[-l,l]
	 * @return		a reference to the element specified by the indices.
	 */
	std::complex<double> & operator() (int r, int m, int l);

	/**
	 * Obtain an element of the expansion.
	 *
	 * @param[in] r		remainder index per l and m quantum number
	 * @param[in] l		angular quantum number l of the spherical harmonic
	 * @param[in] m		magnetic quantum number m of the spherical harmonic, range is m=[-l,l]
	 * @return		a copy of the element specified by the indices.
	 */
	std::complex<double> operator() (int r, int m, int l) const;

	/**
	 * Transform the data in this object using a symmetry operation
	 *
	 * @param[in] sop	The symmertry operation representation.
	 */
	void transform(symmetry::SymmetryOperation const & sop);

	/**
	 * Obtain the number of radial points
	 *
	 * @return	the number of radial points
	 */
	int get_num_radial() const;

	/**
	 * Given a set of arbitrary points, interpolate the data by multiplying the spherical harmonic basis.
	 *
	 * @param[in] coordinates			A list of 3N coordinates, x1,y1,z1,x2,...,zN of N points where the function should be interpolated to.
	 * @param[out] interpolated_data	Resized to fit the N data values of the function at the location of the N points.
	 */
	void interpolate(
			std::vector<double> const & coordinates,
			Auxillary::alignedvector::ZV & interpolated_data) const;

	/**
	 * Given a set of arbitrary points, interpolate the data by multiplying the spherical harmonic basis.
	 *
	 * @param[in] coordinates			A list of 3N coordinates, x1,y1,z1,x2,...,zN of N points where the function should be interpolated to.
	 * @param[out] interpolated_data	Must be allocated to the size N. After calling it contains the N data values of the function at the location of the N points.
	 */
	void interpolate(
			std::vector<double> const & coordinates,
			std::complex<double> * interpolated_data) const;

	/**
	 * Set the internal data from a interpolation produced by the class F.
	 *
	 * We want to represent a function \f$ F({\bf r})\f$ in this object as
	 * \f{eqnarray*}{
	 * F({\bf r}) & = & \sum_{l=0}^{l_{{\rm Max}}}\sum_{m=-l}^{l}f_{lm}(\vert{\bf r}-{\bf r}_{0}\vert)
	 * 					 Y_{l}^{m}(\frac{{\bf r}-{\bf r}_{0}}{\vert{\bf r}-{\bf r}_{0}\vert})
	 * \f}
	 * This method will numerically perform the angular integration to obtain the expansion coefficients:
	 * \f{eqnarray*}{
	 * f_{lm}(r) & = & \int_{0}^{2\pi}{\rm d}\phi\int_{0}^{\pi}{\rm d}\theta{\rm cos}(\theta)
	 * 					F[{\bf r}(r,\theta,\phi)-{\bf r}_{0}]{Y_{l}^{m}}^{\ast}(\theta,\phi)
	 * \f}
	 * at each r in the radial grid. For each of these integrals it uses Algorithms::SphereIntegrator.
	 * The algorithm will pick a rule that integrates 2*LMax exactly. The two is because we expect
	 * data to be of type of a product of two spherical harmonics. If each can have power LMax, we have 2LMax
	 * in the max exponent. Note this is heuristic since the radial data can be anything ...
	 *
	 * @tparam F	A functor that implements the function interpolate() with the same signature and meaning
	 * 				as the one above.
	 *
	 * @param[in] f			A const reference to the object providing the data
	 * @param[in] lmax		The maximal angular momentum we want to represent.
	 * @param[in] rgrid		The radial grid.
	 */
	template<class F>
	void fit_to_data(
			F const & f,
			int lmax,
			RadialGrid rgrid );
private:

	int lmax_ = 0;

	int rMax_ = 0;

	Auxillary::alignedvector::ZV data_;

	RadialGrid rgrid_;

	void set_ylm_conjg_data(
			int l, int m,
			Auxillary::alignedvector::DV const & thetas,
			Auxillary::alignedvector::DV const & phis,
			Auxillary::alignedvector::ZV & ylmData) const;
};

} /* namespace AtomicSite */
} /* namespace elephon */

#include "AtomicSite/SphericalHarmonicExpansion.hpp"
#endif /* ELEPHON_ATOMICSITE_SPHERICALHARMONICEXPANSION_H_ */
