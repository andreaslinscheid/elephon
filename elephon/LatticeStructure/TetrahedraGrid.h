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
#include "LatticeStructure/ExtendedSymmetricGrid.h"
#include <vector>
#include <memory>

namespace elephon
{
namespace LatticeStructure
{
// forward declare
namespace detail{ class Tetrahedra;};

class TetrahedraGrid
{
public:
	typedef detail::Tetrahedra Tetrahedra;

	void initialize(std::shared_ptr<const RegularSymmetricGrid> grid);

	void compute_isosurface_integral_weights() const;

	void get_t_weights( ) const;

	int get_n_tetra() const;

	int get_n_reducible_tetra() const;

	std::vector<Tetrahedra> const get_tetra_list() const;

	std::vector<Tetrahedra> const get_reducible_tetra_list() const;

	std::shared_ptr<const RegularSymmetricGrid> get_grid() const;

	int get_reducible_to_irreducible(int ired) const;

	std::vector<int> const & get_irreducible_to_reducible(int iirred) const;

private:

	std::shared_ptr<const RegularSymmetricGrid> grid_;

	std::vector<Tetrahedra> tetras_;

	std::vector<Tetrahedra> reducibleTetras_;

	std::vector<std::int32_t> reducibleToIrreducible_;

	std::vector<std::vector<std::int32_t>> irreducibleToReducible_;

	void split_cube_insert_tetra(
			RegularSymmetricGrid::GridCube const & cube,
			std::shared_ptr<const ExtendedSymmetricGrid> extendedGrid,
			bool diagonal1,
			std::vector<Tetrahedra> & reducibleTetra,
			std::map<Tetrahedra,std::vector<std::int32_t>> & tetraSet);
};

namespace detail
{
class Tetrahedra
{
public:
	Tetrahedra(
			std::vector<int> cornerIndicesData,
			std::vector<int> cornerIndicesExtended,
			std::vector<int> cornerIndicesExtendedReducible,
			std::shared_ptr<const ExtendedSymmetricGrid> extendedGrid);

	int get_multiplicity() const;

	void set_multiplicity(int m);

	std::vector<int> const & get_corner_indices() const;

	void compute_corner_vectors(
			std::vector<double> & p0,
			std::vector<double> & v123 ) const;

	void compute_corner_points(
			std::vector<double> & p0123 ) const;

	friend bool operator< (Tetrahedra const & t1, Tetrahedra const & t2);

private:

	int multiplicity_ = 1;

	std::vector<int> cornerIndicesData_;

	std::vector<int> cornerIndicesExtended_;

	std::vector<int> cornerIndicesExtendedReducible_;

	std::shared_ptr<const ExtendedSymmetricGrid> extendedGrid_;
};

} /* namespace detail */

} /* namespace LatticeStructure */
} /* namespace elephon */

#endif /* ELEPHON_LATTICESTRUCTURE_TETRAHEDRAGRID_H_ */
