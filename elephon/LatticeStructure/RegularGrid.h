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
#include <set>

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
			LatticeModule lattice);

	int get_np_red() const;

	int get_np_irred() const;

	void get_list_lattice_point_indices(
			std::vector<double> const & points,
			std::vector<int> & irreducibleIndices) const;

	void convert_reducible_irreducible(
			std::vector<int> const& reducibleIndices,
			std::vector<int> & irreducibleIndices) const;

	void convert_reducible_irreducible(
			std::vector<int> const& reducibleIndices,
			std::vector<int> & irreducibleIndices,
			std::vector<int> & symmetryIndexRedToIrred) const;

	int get_xyz_to_reducible(std::vector<int> const & xyzTouple) const;

	bool is_reci() const;
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

	std::vector< std::vector<int> > irredToRredMap_;

	std::vector< std::vector<int> > symIrredToRedMap_;

	void construct_grid_vector_set(
			std::set< RegularGrid::GridPoint > & reducibleSet) const;

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
