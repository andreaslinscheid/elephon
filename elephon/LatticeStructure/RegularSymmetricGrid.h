/*	This file RegularGrid.h is part of elephon.
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
 *  Created on: May 20, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_LATTICESTRUCTURE_REGULARSYMMETRICGRID_H_
#define ELEPHON_LATTICESTRUCTURE_REGULARSYMMETRICGRID_H_

#include "LatticeStructure/Symmetry.h"
#include "LatticeStructure/LatticeModule.h"
#include "LatticeStructure/RegularBareGrid.h"
#include <vector>
#include <map>

namespace elephon
{
namespace LatticeStructure
{

class RegularSymmetricGrid : private RegularBareGrid
{
public:
	typedef RegularBareGrid::GridCube GridCube;
	typedef RegularBareGrid::GridPoint GridPoint;

	using RegularBareGrid::get_vector_direct;
	using RegularBareGrid::get_grid_prec;
	using RegularBareGrid::get_xyz_to_reducible;
	using RegularBareGrid::get_reducible_to_xyz;
	using RegularBareGrid::get_grid_dim;
	using RegularBareGrid::get_grid_shift;
	using RegularBareGrid::compute_grid_cubes_surrounding_nongrid_points;
	using RegularBareGrid::get_lattice;
	using RegularBareGrid::get_list_reducible_lattice_point_indices;
	using RegularBareGrid::interpret_fft_dim_input;

	RegularBareGrid view_bare_grid() const;

	void initialize(
			std::vector<int> dim,
			double gridPrec = 1e-6,
			std::vector<double> shift = {0.0,0.0,0.0},
			Symmetry symmetry = Symmetry(),
			LatticeModule lattice = LatticeModule(),
			std::vector<double> const & constraintIrreducible = std::vector<double>());

	int get_np_red() const;

	int get_np_irred() const;

	void get_list_lattice_point_indices(
			std::vector<double> const & points,
			std::vector<int> & irreducibleIndices) const;

	void convert_reducible_irreducible(
			std::vector<int> const& reducibleIndices,
			std::vector<int> & irreducibleIndices) const;

	std::vector< std::vector<int> > const &  get_maps_irreducible_to_reducible() const;

	std::vector< std::vector<int> > const &  get_maps_sym_irred_to_reducible() const;

	std::vector<int> const &  get_maps_sym_red_to_irreducible() const;

	std::vector<int> const &  get_maps_red_to_irreducible() const;

	bool is_reci() const;

	LatticeStructure::Symmetry const & get_symmetry() const;

private:

	int nIrredPts_ = 0;

	LatticeStructure::Symmetry symmetry_;

	std::vector<int> redToIrredMap_;

	std::vector<int> symRedToIrredMap_;

	std::vector< std::vector<int> > irredToRedMap_;

	std::vector< std::vector<int> > symIrredToRedMap_;

};

} /* namespace LatticeStructure */
} /* namespace elephon */

#endif /* ELEPHON_LATTICESTRUCTURE_REGULARSYMMETRICGRID_H_ */
