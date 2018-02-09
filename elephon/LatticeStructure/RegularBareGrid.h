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
#include "symmetry/SymmetryOperation.h"
#include <vector>
#include <map>

namespace elephon
{
namespace LatticeStructure
{

namespace detail
{
class GridPt;
struct GridCubeImpl;
}; /*namespace detail */

/**
 * Manages a regular grid of 3D vectors in the range for x,y and z of [-0.5,0.0[ in internal coordinates.
 *
 * Internal coordinates are fractions of the lattice basis vectors. See LatticeStructure::LatticeModule.
 * A grid vector has a comparison in the form of a operator< which is done element wise. Any x and x' are
 * said to be qual if they differ by no more than a gridPrec_ value. Then, for two vectors if x is <equal>
 * we compare y and so on.
 */
class RegularBareGrid
{
public:

	typedef struct detail::GridCubeImpl GridCube;

	typedef class detail::GridPt GridPoint;

	/**
	 * Empty constructor, not in a legal state - call initialize() before usage.
	 */
	RegularBareGrid();

	/**
	 *	Construct and directly call initialize()
	 */
	RegularBareGrid(
			std::vector<int> dim,
			bool isReciprocal = false,
			double gridPrec = 1e-6,
			std::vector<double> shift = {0.0,0.0,0.0},
			LatticeModule lattice = LatticeModule());

	/**
	 *	Set the grid.
	 *
	 * @param[in] dim			a 3D vector with the numer of points in each x, y and z direction.
	 * @param[in] isReciprocal	Flag if the grid is sampling reciprocal space.
	 * @param[in] gridPrec		max difference below when to consider two vector components equal
	 * @param[in] shift			a offset of each vector in the grid from the origin. Enter in units
	 * 							of the lattice basis and with a scale of the number of points in each
	 * 							direction. Thus, 0.5, 0, 0.25 will shift every x coordinate of a grid vector by half a lattice spacing
	 * 							in x, does not shift y and shift z by 1/4 of a lattice spacing in z direction.
	 * @param[in] lattice		The lattice module defining the basis.
	 */
	void initialize(
			std::vector<int> dim,
			bool isReciprocal = false,
			double gridPrec = 1e-6,
			std::vector<double> shift = {0.0,0.0,0.0},
			LatticeModule lattice = LatticeModule());

	int get_num_points() const;

	double get_grid_prec() const;

	std::vector<double> get_vector_direct(int i) const;

	/**
	 * Convert a reducible grid index into a grid vector.
	 *
	 * This is the variant of the function that does not reallocate memory
	 * which should be faster where that matters.
	 *
	 * @param i		The input reducible grid index.
	 * @param xyz	a 3 component vector used as buffer. Overwritten on output.
	 * @param v		a 3 component double vector. On output set to the grid reducible vector.
	 */
	void get_vector_direct(int i, std::vector<int> & xyz, std::vector<double> & v) const;

	std::vector<double> get_vectors_direct(std::vector<int> const & reducibleIndices) const;

	std::vector<double> get_all_vectors_grid() const;

	int get_xyz_to_reducible(std::vector<int> const & xyzTouple) const;

	/**
	 * Compute the reducible grid index for given xyz tuple while xyz can have indices in neighboring cells.
	 *
	 * NOTE: only one iteration of periodicty is accounted for. (!!!)
	 *
	 * @param xyzTouple		x y and z coordinates of the grid.
	 * @return				the reducible grid index.
	 */
	int get_xyz_to_reducible_periodic(std::vector<int> xyzTouple) const;

	std::vector<int> get_reducible_to_xyz(int i) const;

	void get_reducible_to_xyz(int i, std::vector<int> & xyz) const;

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

	void transform( symmetry::SymmetryOperation const & sop );

	std::vector<double> get_coords() const;

	friend bool operator< (GridPt const& d1, GridPt const& d2);
private:

	double gridPrec_ = 0;

	std::vector<double> coords_;

	void map_back();
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
