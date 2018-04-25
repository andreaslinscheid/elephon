/*	This file PrimitiveToSupercellConnection.h is part of elephon.
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
 *  Created on: Feb 2, 2018
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_LATTICESTRUCTURE_PRIMITIVETOSUPERCELLCONNECTION_H_
#define ELEPHON_LATTICESTRUCTURE_PRIMITIVETOSUPERCELLCONNECTION_H_

#include "Auxillary/AlignedVector.h"
#include "LatticeStructure/Atom.h"
#include <memory>
#include <vector>
#include <map>

namespace elephon
{
namespace LatticeStructure
{

class UnitCell;
class RegularBareGrid;

/**
 * Handles the transformation of 3D vectors in coordinates of either primitive cell or supercell to the respective other.
 *
 * There are essentially 4 potential cells we worry about.
 *	@todo Input cell must be primitive right now, change that!
 *	There is an input cell (1.) that can be a supercell. We provide functionality (not yet, but later ...) to identify a primitive (2.)
 *	cell in this input cell.
 *	Then, we have a supercell for the vibrational calculation (3.). Furthermore, we have an embedding supercell (4.).
 *	This is the normal supercell for vibrational calculation, (3.), extended to be symmetetrically arranged, so that the primitve cell
 *	is in the middle.  For a regular grid in the primitive cell, the regular grid of the embedding cell is characterized by the property:
 *		1. Each point in Cartesian coordinates corresponds to a primitive grid vector by adding a lattice vector.
 *			In particular for \f$ \bf{R} = \bf{0}\f$ the primitive points are a subset of the embedding grid.
 *		2. Each point has an inverse w.r.p.t. the center which coincides with the primitive cell.
 */
class PrimitiveToSupercellConnection
{
public:

	/**
	 *	Set the internal data.
	 *	Please note: Currently, we only allow transformation where the unit cell is a 3x3x3 supercell.
	 *
	 *	@todo	Allow discovery of more complex relationships between primitive and composite cells.
	 *
	 * @param[in] primitiveCell		A ptr to the primitive cell data
	 * @param[in] superCell			A ptr to the supercell. Must have an integer multiple of primitive volumes, obviously.
	 */
	void initialize(
			std::shared_ptr<const UnitCell> primitiveCell,
			std::shared_ptr<const UnitCell> superCell );

	/**
	 * For an atom index in the supercell, return the atom index of the primitive cell.
	 *
	 * Mathematically, we find \f$ \tau_{\kappa}^{\rm primitive} = \tau_{\kappa^{\prime}}^{\rm SC} + {\bf R}\f$ for suitable \f${\bf R}\f$
	 *
	 * @param iASC	The index of the atom in the supercell.
	 * @return	The primitive atom index mapped to Atom[iASC] by a lattice translation.
	 */
	int get_equiv_atom_primitive(int iASC) const;

	/**
	 * For an atom index in the supercell, return the lattice vector index taking one from the primitive cell to the atom.
	 *
	 * Mathematically, we find \f$ \tau_{\kappa}^{\rm primitive} = \tau_{\kappa^{\prime}}^{\rm SC} + {\bf R}\f$ for suitable \f${\bf R}\f$
	 *
	 * @param[in] iASC	The index of the atom in the supercell.
	 * @return	index of \f$ {\bf R}\f$ with x,y and z coordinate such that
	 * 			\f$ \tau_{\kappa}^{\rm primitive} = \tau_{\kappa^{\prime}}^{\rm SC} + {\bf R}\f$.
	 */
	int get_equiv_atom_primitive_lattice_vector_index(int iASC) const;

	/**
	 * Obtain the atom index of an atom.
	 *
	 * This method returns -1 if the atom does not belong to the lattice.
	 *
	 * @param[in] a					The atom to be found.
	 * @param[in] primitiveCell		If true the atom will be searched in the primitive cell, otherwise we use the supercell.
	 * @return						The index of the atom or -1 if its not found.
	 */
	int find_atom(Atom const & a, bool primitiveCell) const;

	/**
	 * Get the list of supercell vectors and index map.
	 * @param[out] Rvectors			Resizes to the N x 3 supercell vectors.
	 * @param[out] RVectorIndexMap	A multiarray [MaxRVectorX]...[MaxRVectorZ] that returns the index of a vector or -1 if that vector
	 * 								is not in the set.
	 */
	void get_supercell_vectors(Auxillary::Multi_array<int,2> & Rvectors,
			Auxillary::Multi_array<int,3> & RVectorIndexMap) const;

	/**
	 * Get a supercell vector.
	 * @param[in] iR	Index of the vector in the list of supercell vectors.
	 * @return	an array with the 3 components x,y and z of a supercell vector.
	 */
	std::array<int,3> get_supercell_vector(int iR) const;

	/**
	 * Transform a vector represented in the supercell to coordinates of the primitive cell.
	 *
	 * @param[in,out] vec	A 3 component vector that will overwritten with the coordintes in units of the primitive cell.
	 */
	void supercell_to_primitive_coordinates(std::vector<double> & vec) const;

	/**
	 * Transform a vector represented in coorinates of the primitive cell to coordinates of the supercell.
	 * @param[in,out] vec	A 3 component vector that will overwritten with the coordintes in units of the supercell.
	 */
	void primitive_to_supercell_coordinates(std::vector<double> & vec) const;

	/**
	 * Convert an atom index in the primitve cell to a supercell atom index of the equivalent atom.
	 *
	 * 'Equivalent' means that given the definition of the primitive cell w.r.p.t. the supercell
	 * an equivalent atom has the same cartesian coordinates.
	 *
	 * @param primitiveAtomIndex	Integer index referencing the number of the atom in the primitive cell atom list.
	 * @return	The index of the equivalent atom in the list of atoms in the supercell.
	 */
	int primitive_to_supercell_atom_index(int primitiveAtomIndex) const;

	/**
	 * Get size of the supercell
	 * @return	the number of times how often the primitive cell fits into the supercell.
	 */
	int get_supercell_volume_factor() const;

	/**
	 * Generate the data for the embedding supercell.
	 *
	 * See PrimitiveToSupercellConnection documentation for what the embedding cell is.
	 *
	 * @todo: The algorithm to figure out the embedding cell is not optimal. It can be improved by only considering
	 * the primitive cell vectors and the grid dimension. Current implemention loops all grid points which can get heavy...
	 *
	 * @param[in] primitiveCellGrid		Regular grid sampling the primitive cell
	 * @param[in] supercellGrid			Regular grid sampling the supercell
	 * @param[out] minMaxEachDim		Set to the -min (which equals max) of the box in lattice vector space enclosing the embedding cell.
	 * @param[out] listAllRVectors		resized to hold [N_em][3] embedding lattice vector indices and the coordinates.
	 * @param[out] RVectorIndexMap		resized to the 'Box' [-minMaxEachDim[0]:minMaxEachDim[0]][...][-minMaxEachDim[2]:minMaxEachDim[2]] and returns
	 * 									the index of the respective vector index in \p listAllRVectors or -1 if that vector is not in the list
	 * @param[out] embeddingSuperCellGrid
	 */
	void build_embedding_supercell(
			LatticeStructure::RegularBareGrid const & primitiveCellGrid,
			LatticeStructure::RegularBareGrid const & supercellGrid,
			std::array<int,3> & minMaxEachDim,
			Auxillary::Multi_array<int,2> & listRVectors,
			Auxillary::Multi_array<int,3> & RVectorIndexMap,
			LatticeStructure::RegularBareGrid & embeddingSuperCellGrid ) const;

	/**
	 * Find a map from a point in the supercell to the primitive cell centered around a given atom.
	 *
	 * In this method, we consider essentially 3 cells.
	 * 1. The particular primitive cell where a given atom is in the center.
	 * 2. The supercell as computed by the electronic structure code. The grid on which data is represented
	 * 		is shifted such that the given atom is in the center.
	 * 3. The embedding supercell. This is a supercell described by the set of lattice vectors \f${\bf R}\f$ such that
	 * 		\f$ ({\bf r}_{sc}-{\tau}_{0}) = {\bf r}_{0} + {\bf R} \f$
	 * 		Note: If the supercell dimension is not even the grid coordinates of the supercell cut a primitive cell in two.
	 * 		Then, we have vectors \f$ {\bf R} \f$ in units of the shifted primitive cell that lay outside of the supercell.
	 * 		If the supercell is a simple even multiple of the primitve cell in each direction the embedding cell and the
	 * 		supercell are identical.
	 * As a result of this method, we create a table which maps a point index in real space in the supercell
	 * into a point index in the primitive cell plus a lattice vector
	 * For example if there is one atom at (0.25, 0.25, 0.25) in the originial primitive cell, and we have
	 * a supercell of (2,2,2), this shift will be applied
	 * to the entire supercell before the mapping to the primitive cell in the range [-0.5, 0.5[ in each coordinate
	 * will be computed.
	 *	@todo	Allow discovery of more complex relationships between primitive and composite cells.
	 *
	 * @param[in] primitiveCellGrid		The grid which describes the real space of the primitive cell.
	 * 									Taken to sample the primitive cell shifted to every of the atoms in \p atomsUC.
	 * @param[in] supercellGrid			The grid which describes the real space of the supercell. Taken to sample the supercell
	 * 									shifted to every of the atoms in \p atomsUC.
	 * 									The grid which describes the real space of the supercell. Taken to sample the supercell
	 * @param[in] atomsUC				A list of atoms in the primitive cell. An index map and lattice vector map
	 * 									will be computed for each of these atoms, assuming any one of them in the center
	 * 									of a variable primitive cell.
	 * @param[in] minMaxEachDim			The Box enclosing the entire embedding cell.
	 * @param[out] indexMap				a multiarray that will be resized to hold 1) for each atom 2) a list of real space indices
	 * 									referencing the consecutively ordered real space index in the primitive cell
	 * @param[out] latticeVectorMap		a multiarray that will be resized to hold 1) for each atom 2) a list
	 * 									real space index in the supercell and the 3 dimension is the x,y and z coordinate of the R vectors
	 */
	void build_supercell_to_primitive(
			LatticeStructure::RegularBareGrid const & primitiveCellGrid,
			LatticeStructure::RegularBareGrid const & supercellGrid,
			std::vector<LatticeStructure::Atom> const & atomsUC,
			std::array<int,3> const & minMaxEachDim,
			Auxillary::Multi_array<int,2> & indexMap,
			Auxillary::Multi_array<int,3> & latticeVectorMap ) const;

	/**
	 * Get the matrix that transforms vector in the supercell to coordinates of the primitive cell when multiplied from the left.
	 *
	 * @return	A constant reference to the internally stored transformation matrix
	 */
	Auxillary::Multi_array<int,2> const & get_supercell_matrix() const;
private:

	/// When multiplied from the left it transforms a vector in the supercell into coordinates
	/// of the primitive cell (which are of cause not necessarily in the range [-0.5,0.5[)
	Auxillary::Multi_array<int,2> supercellToPrimitiveCoordsMatrix_
										= Auxillary::Multi_array<int,2>(boost::extents[3][3]);

	/// When multiplied from the left it transforms a vector in the primitive cell into coordinates
	/// of the supercell cell
	Auxillary::Multi_array<double,2> primitiveToSupercellCoordsMatrix_
										= Auxillary::Multi_array<double,2>(boost::extents[3][3]);

	/// N x 3 array of the explicit N lattice vectors of the supercell from the POV of the primitive cell.
	Auxillary::Multi_array<int,2> Rvectors_;

	/// An index array returning for given R vector the index if has in the list of R vectors.
	Auxillary::Multi_array<int,3> indexRVectorMap_;

	std::map<LatticeStructure::Atom,int> lookupPrimitive_;

	std::map<LatticeStructure::Atom,int> lookupSupercell_;

	std::shared_ptr<const UnitCell> primitiveCell_;

	std::shared_ptr<const UnitCell> superCell_;

	std::vector<int> superCellToPrimitveAtomIndex_;

	Auxillary::Multi_array<int,2> superCellToPrimitveRVector_;

	void discover_primitive_cell(
			std::shared_ptr<const UnitCell> superCell,
			UnitCell & primitiveCell) const;

	void supercell_to_primitive_coordinates_no_shift(std::vector<double> & vec) const;

	void primitive_to_supercell_coordinates_no_shift(std::vector<double> & vec) const;

	void compute_get_equiv_atom_primitive(int iASC, std::array<int,3> & R, int & iAPC) const;
};

} /* namespace LatticeStructure */
} /* namespace elephon */

#endif /* ELEPHON_LATTICESTRUCTURE_PRIMITIVETOSUPERCELLCONNECTION_H_ */
