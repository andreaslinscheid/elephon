/*	This file TetrahedraGrid.h is part of elephon.
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
 *  Created on: Oct 8, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_LATTICESTRUCTURE_TETRAHEDRAGRID_H_
#define ELEPHON_LATTICESTRUCTURE_TETRAHEDRAGRID_H_

#include "LatticeStructure/RegularSymmetricGrid.h"
#include "LatticeStructure/Tetrahedron.h"
#include <vector>
#include <memory>

namespace elephon
{
namespace LatticeStructure
{

class TetrahedraGrid
{
public:

	void initialize(std::shared_ptr<const RegularSymmetricGrid> grid);

	void compute_isosurface_integral_weights() const;

	void get_t_weights( ) const;

	int get_n_tetra() const;

	int get_n_reducible_tetra() const;

	std::vector<Tetrahedron> const get_tetra_list() const;

	std::vector<Tetrahedron> const get_reducible_tetra_list() const;

	std::shared_ptr<const RegularSymmetricGrid> get_grid() const;

	int get_reducible_to_irreducible(int ired) const;

	std::vector<int> const & get_irreducible_to_reducible(int iirred) const;

	void compute_grid_tetrahedra_surrounding_nongrid_points(
			std::vector<double> const & nonGridPoints,
			std::map<Tetrahedron,std::vector<int>> & tetras) const;
private:

	bool mainDiagon_0_6_ = false;

	std::shared_ptr<const RegularSymmetricGrid> grid_;

	std::shared_ptr<ExtendedSymmetricGrid> extendedGrid_;

	std::vector<Tetrahedron> tetras_;

	std::vector<Tetrahedron> reducibleTetras_;

	std::vector<std::int32_t> reducibleToIrreducible_;

	std::vector<std::vector<std::int32_t>> irreducibleToReducible_;

	void split_cube_insert_tetra(
			RegularSymmetricGrid::GridCube const & cube,
			std::vector<Tetrahedron> & reducibleTetra,
			std::map<Tetrahedron,std::vector<std::int32_t>> & tetraSet) const;
};

} /* namespace LatticeStructure */
} /* namespace elephon */

#endif /* ELEPHON_LATTICESTRUCTURE_TETRAHEDRAGRID_H_ */
