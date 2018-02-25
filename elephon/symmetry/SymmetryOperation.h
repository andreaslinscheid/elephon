/*	This file SymmetryOperation.h is part of elephon.
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
 *  Created on: Jan 16, 2018
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_SYMMETRY_SYMMETRYOPERATION_H_
#define ELEPHON_SYMMETRY_SYMMETRYOPERATION_H_

#include "AtomicSite/ASSymmetry.h"
#include <vector>
#include <memory>

namespace elephon
{

namespace LatticeStructure { class RegularBareGrid; };

namespace symmetry
{

/**
 * Implements functionality of an elementary symmetry operation.
 *
 * The logic for application of a symmetry operation is that this class implements the elementary transformation
 * on a vector or a matrix, while composed objects such as AtomicSite::SphericalHarmonicExpansion implement the
 * full transformation and call the functionality in this object for each component.
 *
 */
class SymmetryOperation
{
public:

	/**
	 * Create the Symmetry operation object from plain data.
	 *
	 * @param[in] pointGroupLatticeBasisBegin	Iterator to the begin of integer point group matrix in the basis of lattice.
	 * @param[in] pointGroupLatticeBasisEnd		Iterator past the last element of the point group matrix in the basis of lattice.
	 * @param[in] fractLatticeBasisBegin		Iterator to the begin of the fractional translation vector in the basis of lattice.
	 * @param[in] fractLatticeBasisEnd			Iterator past the last element of the fractional translation vector in the basis of lattice.
	 * @param[in] pointGroupCarthBasisBegin		Iterator to the begin of point group matrix in Cartesian coordinates.
	 * @param[in] pointGroupCarthBasisEnd		Iterator past the last element of the group matrix in Cartesian coordinates.
	 * @param[in] fractCarthBasisBegin			Iterator to the begin of the fractional translation vector in Cartesian coordinates.
	 * @param[in] fractCarthBasisEnd			Iterator past the last element of the fractional translation vector in Cartesian coordinates.
	 * @param[in] radSym_ptr					shared pointer with the symmetry operation in the basis of spherical harmonics
	 */
	SymmetryOperation(
			std::vector<int>::const_iterator pointGroupLatticeBasisBegin,
			std::vector<int>::const_iterator pointGroupLatticeBasisEnd,
			std::vector<double>::const_iterator fractLatticeBasisBegin,
			std::vector<double>::const_iterator fractLatticeBasisEnd,
			std::vector<double>::const_iterator pointGroupCarthBasisBegin,
			std::vector<double>::const_iterator pointGroupCarthBasisEnd,
			std::vector<double>::const_iterator fractCarthBasisBegin,
			std::vector<double>::const_iterator fractCarthBasisEnd,
			std::shared_ptr<const AtomicSite::ASSymmetry::RadSym> radSym_ptr);

	void rotate( std::vector<int> & v ) const;

	void apply( std::vector<double> & v, bool latticePeriodic = true) const;

	std::shared_ptr<const AtomicSite::ASSymmetry::RadSym>
	get_radial_symmetry_operator() const;

	template<class VT>
	void rotate_cart( VT & v ) const;

	template<typename iterator>
	void rotate_cart( iterator begin, iterator end ) const;

	template<class VT>
	void rotate_matrix_cart(VT & m) const;

	double get_carth_rot_matrix (int i, int j) const;

	double get_carth_frac_trans(int i) const;

	int get_lat_rot_matrix (int i, int j) const;

	int get_lat_frac_trans(int i) const;

	template<class VT, class Grid>
	void transform_field_regular_grid(
			Grid const & grid,
			VT & functionData) const;

	/**
	 * Apply the inverse of the rotation operator to expansion coefficients of spherical harmonics.
	 *
	 *	Since the symmetry operation acts on the function, not the coefficients, after a transformation,
	 *	the coefficient belongs to target of the transformation. The relevant derivation is
	 *	 \f{eqnarray*}{
	 *		 A(\bf{r}) & = & \sum_{l=0}^{l_{{\rm Max}}}\sum_{m=-l}^{l}a_{l,m}(\vert{\bf{r}}\vert)Y_{l}^{m}(\frac{\bf{r}}{\vert\bf{r}\vert}) \\
	 *		 \hat{R}\cdot A(\bf{r}) & = &\sum_{l=0}^{l_{{\rm Max}}}\sum_{m,m^{\prime}=-l}^{l}
	 *		 									a_{l,m}(\vert{\bf{r}}\vert)D_{mm^{\prime}}^{l\ast}(\alpha,\beta,\gamma)Y_{l}^{m^{\prime}}(\frac{{\bf{r}}}{\vert\bf{r}\vert})\\
	 *		 					&& \equiv A^{\prime}({\bf{r}})=\sum_{l=0}^{l_{{\rm Max}}}\sum_{m=-l}^{l}a_{l,m}^{\prime}(\vert{{\bf{r}}}\vert)
	 *		 					              Y_{l}^{m}(\frac{{\bf{r}}}{\vert{\bf{r}}\vert})
	 *	 \f}
	 *	 such that, by comparing linearly independent coefficients, we identify
	 *	 \f{eqnarray*}{
	 *	 	a_{l,m}^{\prime}(\vert{\bf{r}}\vert)=\sum_{m^{\prime}=-l}^{l}a_{l,m^{\prime}}(\vert{\bf r}\vert)D_{m^{\prime}m}^{l\ast}(\alpha,\beta,\gamma)
	 *	  \f}
	 *	  This method computes the array \f$a_{l,m}^{\prime}\f$ from the array \f$a_{l,m}\f$
	 *
	 * @tparam	VT	Intended to be a vector of double, float or their complex variants with custom allocator.
	 *
	 * @param[in] lMax						The maximal l that the expansion includes. This it defines the range l=[0,lMax].
	 * @param[in] numDataPerLM				For each l and m=[=l,l] the data to be transformed as numDataPerLM elements.
	 * @param[in,out] dataToBeTransformed	Data with the expansion coefficients and a layout of l,m according to AtomicSite::SphericalHarmonicExpansion::angular_momentum_layout
	 *										where for each l,m there numDataPerLM data elements to be transformed as a block.
	 */
	template<class VT>
	void rotate_radial_data(
			int lMax,
			int numDataPerLM,
			VT & dataToBeTransformed) const;

	/**
	 * For the documentation, please see rotate_radial_data()
	 */
	template<class interator>
	void rotate_radial_data(
			int lMax,
			int numDataPerLM,
			interator dataToBeTransformedBegin,
			interator dataToBeTransformedEnd) const;

private:

	int ptgroup[9];

	double ptgCart[9];

	double fracTrans[3];

	double fracTransCart[3];

	std::shared_ptr<const AtomicSite::ASSymmetry::RadSym> radSym_ptr_;
};

} /* namespace symmetry */
} /* namespace elephon */

#include "symmetry/SymmetryOperation.hpp"
#endif /* ELEPHON_SYMMETRY_SYMMETRYOPERATION_H_ */
