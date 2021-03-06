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
	 * Set the internal data from expansion data within real spherical harmonics.
	 *
	 * The potential is described e.g. in VASP in terms of real spherical harmonics
	 * \f{eqnarray*}
	 * 	\tilde{Y}_{l}^{m}(\theta,\phi)
	 * 		& = & \begin{cases}
	 *			\frac{{\rm i}}{\sqrt{2}}[Y_{l}^{m}(\theta,\phi)-(-1)^{m}Y_{l}^{-m}(\theta,\phi)] & m<0\\
	 *			Y_{l}^{m}(\theta,\phi) & m=0\\
	 *			\frac{1}{\sqrt{2}}[Y_{l}^{-m}(\theta,\phi)+(-1)^{m}Y_{l}^{m}(\theta,\phi)] & m>0
	 *			\end{cases}
	 * \f}
	 *	 as the series
	 * \f{eqnarray*}{
	 *	 v_{{\rm scf}}^{{\rm AE}-1}({\bf r}) & = & \sum_{lm}\tilde{Y}_{l}^{m}({\bf r}/\vert{\bf r}\vert)v_{lm}^{{\rm R}}(\vert{\bf r}\vert)
	 * \f}
	 *	 which can be written, again, in terms of the complex spherical harmonics as we need in this class via
	 * \f{eqnarray*}{
	 * v_{{\rm scf}}^{{\rm AE}-1}(\boldsymbol{r})
	 * &=&	\sum_{l,m<0}\frac{{\rm i}}{\sqrt{2}}[Y_{l}^{m}(\theta,\phi)-(-1)^{m}Y_{l}^{-m}(\theta,\phi)]v_{lm}^{{\rm R}}(\vert\boldsymbol{r}\vert) \\
	 * &&	+\sum_{l}Y_{l}^{0}(\theta,\phi)v_{l0}^{{\rm R}}(\vert\boldsymbol{r}\vert)\\
	 * &&	+\sum_{l,m>0}\frac{1}{\sqrt{2}}[Y_{l}^{-m}(\theta,\phi)+(-1)^{m}Y_{l}^{m}(\theta,\phi)]v_{lm}^{{\rm R}}(\vert\boldsymbol{r}\vert)\\
	 * &=&	\sum_{l}Y_{l}^{0}(\theta,\phi)v_{l0}^{{\rm R}}(\vert\boldsymbol{r}\vert)+\sum_{l,m>0}\frac{1}{\sqrt{2}}[Y_{l}^{-m}(\theta,\phi){\rm i}
	 * 			v_{l,-m}^{{\rm R}}(\vert\boldsymbol{r}\vert)-(-1)^{-m}Y_{l}^{m}(\theta,\phi){\rm i}v_{l,-m}^{{\rm R}}(\boldsymbol{\vert\boldsymbol{r}}\vert)\\
	 * &&	+Y_{l}^{-m}(\theta,\phi)v_{lm}^{{\rm R}}(\vert\boldsymbol{r}\vert)+(-1)^{m}Y_{l}^{m}(\theta,\phi)v_{lm}^{{\rm R}}(\vert\boldsymbol{r}\vert)]\\
	 * &=&	\sum_{l}Y_{l}^{0}(\theta,\phi)v_{l}^{0}(\vert\boldsymbol{r}\vert)+\sum_{l,m<0}Y_{l}^{m}(\theta,\phi)\frac{1}{\sqrt{2}}
	 * 			[v_{lm}^{{\rm R}}(\vert\boldsymbol{r}\vert)+{\rm i}v_{l,-m}^{{\rm R}}(\vert\boldsymbol{r}\vert)]\\
	 * &&	+\sum_{l,m>0}Y_{l}^{m}(\theta,\phi)\frac{(-1)^{m}}{\sqrt{2}}[v_{lm}^{{\rm R}}(\vert\boldsymbol{r}\vert)
	 * 			-{\rm i}v_{l,-m}^{{\rm R}}(\vert\boldsymbol{r}\vert)]\\
	 * &&	\equiv\sum_{l=0}^{l_{max}}\sum_{m=-l}^{l}Y_{l}^{m}(\theta,\phi)v_{l,m}^{{\rm AE}}(\vert\boldsymbol{r}\vert)	 v_{{\rm scf}}^{{\rm AE}-1}({\bf r})
	 * \f}
	 *	 with
	 * \f{eqnarray*}{
	 * v_{l,m}^{{\rm AE}}(\vert{\bf r}\vert)
	 * 	& = & \begin{cases}
	 *		\frac{1}{\sqrt{2}}[v_{l,-m}^{{\rm R}}(\vert{\bf r}\vert)+{\rm i}v_{lm}^{{\rm R}}(\vert{\bf r}\vert)] & m<0	\\
	 *		v_{l0}^{{\rm R}}(\boldsymbol{\vert{\bf r}}\vert) & m=0	\\
	 *		\frac{(-1)^{m}}{\sqrt{2}}[v_{lm}^{{\rm R}}(\vert{\bf r}\vert)-{\rm i}v_{l,-m}^{{\rm R}}(\vert{\bf r}\vert)] & m>0
	 *	\end{cases}
	 * \f}
	 *
	 *	@tparam VT a container type with linear storage such as vector
	 *
	 * @param[in] lmax					maximal angular momentum that occurs l=(0,...,lmax)
	 * @param[in] dataExpansionRealSH	Data with the expansion data in terms of real spherical harmonics. Layout is the same as in initialize()
	 * @param[in] rgrid					The radial grid.
	 */
	template<class VT>
	void initialize_from_real(
			int lmax,
			VT const & dataExpansionRealSH,
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
	 * Change the center of the expansion.
	 *
	 * @param[in] center	The 3 x,y and z components of the new center.
	 */
	void set_center(std::vector<double> center);

	/**
	 * Obtain the radial grid
	 *
	 * @return	a constant reference to the radial grid.
	 */
	RadialGrid const & get_radial_grid() const;

	/**
	 * Get the maximal contained angular momentum number.
	 * @return	The highest angular momentum number for which we have data.
	 */
	int get_l_max() const;

	/**
	 * Given a set of arbitrary points, interpolate the data by multiplying the spherical harmonic basis.
	 *
	 * @tparam VT container type with linear storage such as vector
	 *
	 * @param[in] coordinates			A list of 3N coordinates, x1,y1,z1,x2,...,zN of N points where the function should be interpolated to.
	 * @param[out] interpolated_data	Resized to fit the N data values of the function at the location of the N points.
	 */
	template<class VT>
	void interpolate(
			VT const & coordinates,
			Auxillary::alignedvector::ZV & interpolated_data) const;

	/**
	 * Given a set of arbitrary points, interpolate the data by multiplying the spherical harmonic basis.
	 *
	 * @tparam VT container type with linear storage such as vector
	 *
	 * @param[in] coordinates			A list of 3N coordinates, x1,y1,z1,x2,...,zN of N points where the function should be interpolated to.
	 * @param[out] interpolated_data	Must be allocated to the size N. After calling it contains the N data values of the function at the location of the N points.
	 */
	template<class VT>
	void interpolate(
			VT const & coordinates,
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

	/**
	 * Random access iterator to the data in this container.
	 * @return	beginning of the range.
	 */
	Auxillary::alignedvector::ZV::iterator begin();

	/**
	 * Random access iterator to the end in this container.
	 * @return	one after the last element in the range.
	 */
	Auxillary::alignedvector::ZV::iterator end();

	/**
	 * Constant random access iterator to the data in this container.
	 * @return	beginning of the range.
	 */
	Auxillary::alignedvector::ZV::const_iterator begin() const;

	/**
	 * Constant random access iterator to the end in this container.
	 * @return	one after the last element in the range.
	 */
	Auxillary::alignedvector::ZV::const_iterator end() const;
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
