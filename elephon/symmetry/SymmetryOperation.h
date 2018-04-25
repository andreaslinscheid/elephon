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
	void transform_scalar_field_regular_grid(
			Grid const & grid,
			VT & functionData) const;

	/**
	 * Apply this symmetry operation to vector field data.
	 *
	 * @param[in] grid					Regular grid on which the function data is defined.
	 * @param[in,out] functionData		If transposeFunctionData == false, then we expect for each grid point 3 consecutive values
	 * 									with the vector x,y and z data. If transposeFunctionData == true, then, instead, we expect
	 * 									for each x, y and z the entire block of data for each grid value.
	 * @param[in] transposeFunctionData	Choose if component or grid index is the fast running dimension in the data layout.
	 */
	template<class VT, class Grid>
	void transform_vector_field_regular_grid(
			Grid const & grid,
			VT & functionData,
			bool transposeFunctionData = false) const;

	/**
	 * Apply the inverse of the rotation operator to expansion coefficients of spherical harmonics.
	 *
	 *	Since the symmetry operation acts on the function, not the coefficients, after a transformation,
	 *	the coefficient belongs to target of the transformation. The relevant derivation is
	 *	 \f{eqnarray*}{
	 *		 A(\bf{r}) & = & \sum_{l=0}^{l_{{\rm Max}}}\sum_{m=-l}^{l}a_{l,m}(\vert{\bf{r}}\vert)Y_{l}^{m}(\frac{\bf{r}}{\vert\bf{r}\vert}) \\
	 *		 \hat{R}\cdot A(\bf{r}) & = &\sum_{l=0}^{l_{{\rm Max}}}\sum_{m,m^{\prime}=-l}^{l}
	 *		 									a_{l,m}(\vert{\bf{r}}\vert)D_{m^{\prime}m}^{l}(\alpha,\beta,\gamma)Y_{l}^{m^{\prime}}(\frac{{\bf{r}}}{\vert\bf{r}\vert})\\
	 *		 					&& \equiv A^{\prime}({\bf{r}})=\sum_{l=0}^{l_{{\rm Max}}}\sum_{m=-l}^{l}a_{l,m}^{\prime}(\vert{{\bf{r}}}\vert)
	 *		 					              Y_{l}^{m}(\frac{{\bf{r}}}{\vert{\bf{r}}\vert})
	 *	 \f}
	 *	 such that, by comparing linearly independent coefficients, we identify
	 *	 \f{eqnarray*}{
	 *	 	a_{l,m}^{\prime}(\vert{\bf{r}}\vert)=\sum_{m^{\prime}=-l}^{l}D_{mm^{\prime}}^{l}(\alpha,\beta,\gamma)a_{l,m^{\prime}}(\vert{\bf r}\vert)
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
	void rotate_radial_scalar_data(
			int lMax,
			int numDataPerLM,
			VT & dataToBeTransformed) const;

	/**
	 * For the documentation, please see rotate_radial_data()
	 */
	template<class interator>
	void rotate_radial_scalar_data(
			int lMax,
			int numDataPerLM,
			interator dataToBeTransformedBegin,
			interator dataToBeTransformedEnd) const;

	/**
	 * Apply the inverse of the rotation operator to expansion coefficients of spherical harmonics.
	 *
	 * See rotate_radial_scalar_data() for details on the derivation on why we take in the inverse symmetry operation
	 * that holds similarly for this vector data.
	 * However, if the function is vector valued, we find
	 *	 \f{eqnarray*}{
	 *	 \hat{R}\cdot{\bf B}({\bf r}) & = & \sum_{l=0}^{l_{{\rm Max}}}\sum_{m=-l}^{l}g\cdot{\bf b}(_{lm}(\vert{\bf r}\vert)
	 *	 									\hat{R}\cdot Y_{l}^{m}(\frac{{\bf r}}{\vert{\bf r}\vert})\,,
	 *	  \f}
 	 * and we find by the same analysis
	 *	 \f{eqnarray*}{
	 *	 {\bf b}_{lm}^{\prime}(\vert{\bf r}\vert) & = & g\cdot\sum_{m^{\prime}=-l}^{l}{\bf b}_{lm^{\prime}}(\vert{\bf r}\vert)
	 *	 												D_{m^{\prime}m}^{l\ast}(\alpha,\beta,\gamma)\,.
	 *	  \f}
	 *
	 * @param[in] lMax							The maximal l that the expansion includes. This it defines the range l=[0,lMax] and m=[-l,l] for each l.
	 * @param[in] numDataPerLM					For each l and m=[=l,l] the data to be transformed as numDataPerLM elements.
	 * @param[in,out] dataToBeTransformedBegin	Data with the expansion coefficients, first for x, then for y and finally for z, each channel
	 * 											with a layout of l,m according to AtomicSite::SphericalHarmonicExpansion::angular_momentum_layout
	 * @param[in] dataToBeTransformedEnd		Iterator pointing to after the end of the data block.
	 */
	template<class interator>
	void rotate_radial_vector_data(
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

	template<class Grid>
	void compute_real_space_map(Grid const & grid,
			std::vector<int> & inverseSymOpMap) const;
};

} /* namespace symmetry */
} /* namespace elephon */

#include "symmetry/SymmetryOperation.hpp"
#endif /* ELEPHON_SYMMETRY_SYMMETRYOPERATION_H_ */
