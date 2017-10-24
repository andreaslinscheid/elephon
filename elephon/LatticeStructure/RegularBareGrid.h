/*	This file RegularBareGrid.h is part of elephon.
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
 *  Created on: Jul 7, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_LATTICESTRUCTURE_REGULARBAREGRID_H_
#define ELEPHON_LATTICESTRUCTURE_REGULARBAREGRID_H_

#include "LatticeStructure/LatticeModule.h"
#include "LatticeStructure/Symmetry.h"
#include <vector>
#include <map>

namespace elephon
{
namespace LatticeStructure
{

namespace detail
{
class GridPt;
class GridCubeImpl;
}; /*namespace detail */

class RegularBareGrid
{
public:

	typedef class detail::GridCubeImpl GridCube;

	typedef class detail::GridPt GridPoint;

	void initialize(
			std::vector<int> dim,
			bool isReciprocal = false,
			double gridPrec = 1e-6,
			std::vector<double> shift = {0.0,0.0,0.0},
			LatticeModule lattice = LatticeModule());

	int get_num_points() const;

	double get_grid_prec() const;

	std::vector<double> get_vector_direct(int i) const;

	int get_xyz_to_reducible(std::vector<int> const & xyzTouple) const;

	int get_xyz_to_reducible_periodic(std::vector<int> xyzTouple) const;

	std::vector<int> get_reducible_to_xyz(int i) const;

	std::vector<int> const & get_grid_dim() const;

	std::vector<double> get_grid_shift() const;

	std::vector<int> compute_reducible_cube_indices_surrounding_nongrid_point(
			std::vector<double> const& v) const;

	void find_closest_reducible_grid_points(
			std::vector<double> const & nonGridPoints,
			std::vector<int> & reducibleGridIndices) const;

	void compute_grid_cubes_surrounding_nongrid_points(
			std::vector<double> const & nonGridPoints,
			std::vector<int> & nonGridPtToCubeMap,
			std::vector<GridCube> & cubes) const;

	void get_list_reducible_lattice_point_indices(
			std::vector<double> const & points,
			std::vector<int> & reducibleIndices) const;

	void construct_grid_vector_set(
			std::map< RegularBareGrid::GridPoint, int > & reducibleSet) const;

	LatticeStructure::LatticeModule const & get_lattice() const;

	/**
	 * Since input settings for grid dims require sensible defaults this function handles input interpretation.
	 *
	 * @param fftDim	see documentation of the input option fftd.
	 * @return			A vector of dimension 3 with the interpreted grid dimension.
	 */
	std::vector<int> interpret_fft_dim_input(std::vector<int> fftDim) const;

	void get_grid_cubes(std::vector<GridCube> & cubes) const;
private:

	int numPoints_ = 0;

	bool isReciprocalGrid_ = true;

	double gridPrec_ = 0;

	std::vector<int> pointMesh_{0,0,0};

	std::vector<double> pointShift_{0.0,0.0,0.0};

	LatticeStructure::LatticeModule lattice_;
};

namespace detail
{
class GridPt
{
public:
	GridPt();

	GridPt(std::vector<double> coords, double gridPrec = 1e-6);

	void initialize(std::vector<double> coords, double gridPrec = 1e-6);

	void transform( Symmetry::Sop const & sop );

	friend bool operator< (GridPt const& d1, GridPt const& d2);
private:

	double gridPrec_ = 0;

	std::vector<double> coords_;
};

struct GridCubeImpl
{
	GridCubeImpl( std::vector<int> cornerIndices )
		: cornerIndices_(std::move(cornerIndices)) { };

	std::vector<int> cornerIndices_;

	mutable std::vector<int> containedIrregularPts_;

	friend bool operator< (GridCubeImpl const& d1, GridCubeImpl const& d2);
};

} /* namespace detail */
} /* namespace LatticeStructure */
} /* namespace elephon */

#endif /* ELEPHON_LATTICESTRUCTURE_REGULARBAREGRID_H_ */
