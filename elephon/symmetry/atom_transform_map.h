/*	This file atom_transform_map.h is part of elephon.
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
 *  Created on: Jan 15, 2018
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_SYMMETRY_ATOM_TRANSFORM_MAP_H_
#define ELEPHON_SYMMETRY_ATOM_TRANSFORM_MAP_H_

#include "LatticeStructure/Atom.h"
#include "LatticeStructure/Symmetry.h"
#include <vector>

namespace elephon
{
namespace symmetry
{
/** @file */

/**
 * Compute the atom index change due to the transformation of a symmetry group.
 *
 * @param[in] atoms			A vector with the atoms to be transformed.
 * @param[in] symmetry		The symmetry group of that will be applied.
 * @param[out] rotAtomsMap	Overwritten with a vector of size of the symmetry group, where each element is a vector that
 * 							tells for a given atom index where it was before application of the symmetry operation.
 */
void
atom_transform_map(
		std::vector<LatticeStructure::Atom> const & atoms,
		LatticeStructure::Symmetry const & symmetry,
		std::vector< std::vector<int> > & rotAtomsMap);

/**
 * Compute the atom index change due to the transformation of a symmetry group.
 *
 * @param[in] atoms			A vector with the atoms to be transformed.
 * @param[in] symmetry		The symmetry group of that will be applied.
 * @param[out] rotAtomsMap	Overwritten with a vector of size of the symmetry group, where is element is a vector that
 * 							tells for a given atom index where it is after application of the symmetry operation.
 */
void
atom_transform_map_inverse(
		std::vector<LatticeStructure::Atom> const & atoms,
		LatticeStructure::Symmetry const & symmetry,
		std::vector< std::vector<int> > & rotAtomsMap);

} /* namespace Symmetry */
} /* namespace elephon */

#endif /* ELEPHON_SYMMETRY_ATOM_TRANSFORM_MAP_H_ */
