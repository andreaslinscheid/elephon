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

/**
 * Handles the transformation of 3D vectors in coordinates of either primitive cell or supercell to the respective other.
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
	 * Enter a vector in the cell and get the vector in the primitive cell plus Lattice vector.
	 *
	 * Let the vector in the supercell be \f$ {\bf r}\f$, then \f$ {\bf r} = {\bf r_{\rm UC}} + {\bf R}\f$ .
	 *
	 * @param scVec					A 3 component vector x y and z
	 * @param primVecPlusLatticeVec	A pair with the 3 component vector \f${\bf r_{\rm UC}}\f$ first and
	 * 								the 3 component vector \f${\bf R}\f$ second.
	 */
	void supercell_to_primitive_plus_lattice_vector(
			std::vector<double> const & scVec,
			std::pair<std::vector<double>, std::vector<int>> & primVecPlusLatticeVec) const;

	/**
	 * For an atom index in the supercell, return the atom index of the primitive cell.
	 *
	 * Mathematically, we find \f$ \tau_{\kappa}^{\rm primitive} = \tau_{\kappa^{\prime}}^{\rm SC} + {\bf R}\f$ for suitable \f${\bf R}\f$
	 *
	 * @param iASC	The index of the atom in the supercell.
	 * @return	The primitive atom index mapped to Atom[iASC] by a lattice translation.
	 */
	int equiv_atom_primitive(int iASC) const;

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

	void supercell_vectors(Auxillary::Multi_array<int,2> & Rvectors) const;

	void supercell_to_primitive_coordinates(std::vector<double> & vec) const;

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
	int supercell_volume_factor() const;

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

	std::map<LatticeStructure::Atom,int> lookupPrimitive_;

	std::map<LatticeStructure::Atom,int> lookupSupercell_;

	std::vector<LatticeStructure::Atom> primitiveAtomList_;

	void discover_primitive_cell(
			std::shared_ptr<const UnitCell> superCell,
			UnitCell & primitiveCell) const;
};

} /* namespace LatticeStructure */
} /* namespace elephon */

#endif /* ELEPHON_LATTICESTRUCTURE_PRIMITIVETOSUPERCELLCONNECTION_H_ */
