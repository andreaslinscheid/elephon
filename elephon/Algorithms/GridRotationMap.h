/*	This file GridRotationMap.h is part of elephon.
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
 *  Created on: Sep 26, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_ALGORITHMS_GRIDROTATIONMAP_H_
#define ELEPHON_ALGORITHMS_GRIDROTATIONMAP_H_

#include "LatticeStructure/Symmetry.h"
#include "LatticeStructure/RegularBareGrid.h"
#include <vector>

namespace elephon
{
namespace Algorithms
{

/**
 * Compute an index map for each symmetry operation.
 *
 * @param shift			(input) Displace the regular grid *before* applying the symmetry operation.
 * 						Must be a 3 element vector in grid internal units in the 1st Wigner Seitz cell [-0.5,0.5[
 * 						Must be the coordinate of a grid vectors, such that the grid is mapped to itself within
 * 						periodic boundary conditions.
 * @param grid			(input) The regular grid to be rotated.
 * @param symmetry		(input) The symmetry group who's operations are applied to the grid. Transformations will be handled
 * 						modulo a reciprocal lattice vector.
 * @param rotMap		The transformation in terms of grid indices.
 * 						On output, it will be a vector of vectors of int. For each symmetry operation, it will contain a vector
 * 						of #gridpoint grid indices such that rotMap[isym][ir] will tell from which index 'ir' originated from
 * 						before the symmetry operation was applied.
 */
void compute_grid_rotation_map(
		std::vector<double> const & shift,
		LatticeStructure::RegularBareGrid const & grid,
		LatticeStructure::Symmetry const & symmetry,
		std::vector< std::vector<int> > & rotMap);

} /* namespace Algorithms */
} /* namespace elephon */

#endif /* ELEPHON_ALGORITHMS_GRIDROTATIONMAP_H_ */
