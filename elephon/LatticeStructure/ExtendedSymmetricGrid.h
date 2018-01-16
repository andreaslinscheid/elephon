/*	This file ExtendedSymmetricGrid.h is part of elephon.
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
 *  Created on: Oct 21, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_LATTICESTRUCTURE_EXTENDEDSYMMETRICGRID_H_
#define ELEPHON_LATTICESTRUCTURE_EXTENDEDSYMMETRICGRID_H_

#include "LatticeStructure/RegularSymmetricGrid.h"

namespace elephon
{
namespace LatticeStructure
{

namespace detail{ class GridPoint; };

/**
 * TODO: Problems: The underlying grid point is not strictly in the 1. unit cell. This is
 * 		 different from RegularSymmetricGrid.
 * 		 In order to be able to replace RegularSymmetricGrid::initialize with the version
 * 		 that treats the grid point according to the requirements here, ExtendedSymmetricGrid
 * 		 is declared a friend of RegularSymmetricGrid (see there). It would be good to disentangle
 * 		 the grids in the code by some further abstractification.
 */
class ExtendedSymmetricGrid : private RegularSymmetricGrid
{
public:
	typedef RegularSymmetricGrid::GridCube GridCube;
	typedef detail::GridPoint GridPoint;

	using RegularSymmetricGrid::get_grid_prec;
	using RegularSymmetricGrid::get_xyz_to_reducible;
	using RegularSymmetricGrid::get_reducible_to_xyz;
	using RegularSymmetricGrid::get_xyz_to_reducible_periodic;
	using RegularSymmetricGrid::get_grid_dim;
	using RegularSymmetricGrid::get_grid_shift;
	using RegularSymmetricGrid::get_lattice;
	using RegularSymmetricGrid::get_maps_red_to_irreducible;

	void initialize( 	std::vector<int> dim,
						double gridPrec,
						std::vector<double> shift,
						Symmetry symmetry,
						LatticeModule lattice );

	/**
	 * Compute a 3 vector in units of the (reciprocal) lattice.
	 *
	 * Since this is the extended grid, the largest index in each direction x,y and z
	 * will create a vectors that is not in the first unit cell [0,1[.
	 *
	 * @param reducibleIndex 	consequtively ordered index of the reducible grid.
	 * @return					A vector with 3 components, x, y and z.
	 */
	std::vector<double> get_vector_direct(int reducibleIndex) const;

private:

	std::vector<int> irreducibleGridVectorIndices_;
};

namespace detail
{

/**
 * This grid point is defined in the cell [0,1] in all x, y and z
 * directions. This includes both borders of the cell, the one at 0
 * and the one at 1. Thus, a GridPoint coordinate will be transformed
 * back to the range [0,1] only if it exceeds 1+gridPrec. Values
 * not mapped back but larger (smaller) than 1 (0) will be set to 1 (0).
 */
class GridPoint
{
public:
	GridPoint();

	GridPoint(std::vector<double> coords, std::vector<double> extended, double gridPrec = 1e-6);

	void initialize(std::vector<double> coords, std::vector<double> extended, double gridPrec = 1e-6);

	void transform( symmetry::SymmetryOperation const & sop );

	friend bool operator< (GridPoint const& d1, GridPoint const& d2);

private:

	double gridPrec_ = 0;

	std::vector<double> coords_;

	std::vector<double> extended_;

	void map_back();
};
}; /* namespace detail */

} /* namespace LatticeStructure */
} /* namespace elephon */

#endif /* ELEPHON_LATTICESTRUCTURE_EXTENDEDSYMMETRICGRID_H_ */
