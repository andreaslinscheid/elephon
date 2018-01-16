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

#include <vector>

namespace elephon
{
namespace symmetry
{

/**
 * Implements functionality of an elementary symmetry operation.
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
	 */
	SymmetryOperation(
			std::vector<int>::const_iterator pointGroupLatticeBasisBegin,
			std::vector<int>::const_iterator pointGroupLatticeBasisEnd,
			std::vector<double>::const_iterator fractLatticeBasisBegin,
			std::vector<double>::const_iterator fractLatticeBasisEnd,
			std::vector<double>::const_iterator pointGroupCarthBasisBegin,
			std::vector<double>::const_iterator pointGroupCarthBasisEnd,
			std::vector<double>::const_iterator fractCarthBasisBegin,
			std::vector<double>::const_iterator fractCarthBasisEnd );

	void rotate( std::vector<int> & v ) const;

	void apply( std::vector<double> & v, bool latticePeriodic = true) const;

	void rotate_cart( std::vector<double> & v ) const;

	template<class VT>
	void rotate_matrix_cart(VT & m) const;

	double get_carth_rot_matrix (int i, int j) const;

	double get_carth_frac_trans(int i) const;

	int get_lat_rot_matrix (int i, int j) const;

	int get_lat_frac_trans(int i) const;

private:

	int ptgroup[9];

	double ptgCart[9];

	double fracTrans[3];

	double fracTransCart[3];
};

} /* namespace symmetry */
} /* namespace elephon */

#include "symmetry/SymmetryOperation.hpp"
#endif /* ELEPHON_SYMMETRY_SYMMETRYOPERATION_H_ */
