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

#ifndef ELEPHON_LATTICESTRUCTURE_REGULARGRID_H_
#define ELEPHON_LATTICESTRUCTURE_REGULARGRID_H_

#include "LatticeStructure/Symmetry.h"
#include "LatticeStructure/LatticeModule.h"
#include <vector>
#include <map>

namespace elephon
{
namespace LatticeStructure
{

namespace detail
{
class GridPt;
}; /*namespace detail */

class RegularGrid
{
public:

	typedef class detail::GridPt GridPoint;

	void initialize(
			double gridPrec,
			std::vector<int> dim,
			std::vector<double> shift,
			Symmetry symmetry,
			LatticeModule lattice,
			std::vector<double> const & constraintIrreducible = std::vector<double>());

	int get_np_red() const;

	int get_np_irred() const;

	double get_grid_prec() const;

	std::vector<double> get_vector_direct(int i) const;

	void get_list_lattice_point_indices(
			std::vector<double> const & points,
			std::vector<int> & irreducibleIndices) const;

	void get_list_reducible_lattice_point_indices(
			std::vector<double> const & points,
			std::vector<int> & reducibleIndices) const;

	void convert_reducible_irreducible(
			std::vector<int> const& reducibleIndices,
			std::vector<int> & irreducibleIndices) const;

	std::vector< std::vector<int> > const &  get_maps_irreducible_to_reducible() const;

	std::vector< std::vector<int> > const &  get_maps_sym_irred_to_reducible() const;

	std::vector<int> const &  get_maps_sym_red_to_irreducible() const;

	std::vector<int> const &  get_maps_red_to_irreducible() const;

	int get_xyz_to_reducible(std::vector<int> const & xyzTouple) const;

	std::vector<int> get_reducible_to_xyz(int i) const;

	bool is_reci() const;

	LatticeStructure::Symmetry const & get_symmetry() const;

	std::vector<int> const & get_grid_dim() const;

	std::vector<double> get_grid_shift() const;
private:

	double gridPrec_ = 0;

	int nRedPts_ = 0;

	int nIrredPts_ = 0;

	std::vector<int> pointMesh_;

	std::vector<double> pointShift_;

	LatticeStructure::Symmetry symmetry_;

	LatticeStructure::LatticeModule lattice_;

	std::vector<int> redToIrredMap_;

	std::vector<int> symRedToIrredMap_;

	std::vector< std::vector<int> > irredToRedMap_;

	std::vector< std::vector<int> > symIrredToRedMap_;

	void construct_grid_vector_set(
			std::map< RegularGrid::GridPoint, int > & reducibleSet) const;

};

namespace detail
{
class GridPt
{
public:
	GridPt();

	GridPt(std::vector<double> coords, double gridPrec = 1e-6);

	void initialize(std::vector<double> coords, double gridPrec = 1e-6);

	void transform( LatticeStructure::Symmetry::SymmetryOperation const & op);

	friend bool operator< (GridPt const& d1, GridPt const& d2);
private:

	double gridPrec_ = 0;

	std::vector<double> coords_;
};
}

} /* namespace LatticeStructure */
} /* namespace elephon */

#endif /* ELEPHON_LATTICESTRUCTURE_REGULARGRID_H_ */
