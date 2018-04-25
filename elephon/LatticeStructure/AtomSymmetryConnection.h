/*	This file AtomSymmetryConnection.h is part of elephon.
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
 *  Created on: Feb 12, 2018
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_LATTICESTRUCTURE_ATOMSYMMETRYCONNECTION_H_
#define ELEPHON_LATTICESTRUCTURE_ATOMSYMMETRYCONNECTION_H_

#include <vector>

namespace elephon
{
namespace LatticeStructure
{

class Atom;
class Symmetry;

/**
 * Stores symmetry information about atoms in a unit cell.
 */
class AtomSymmetryConnection
{
public:

	/**
	 * Initialize the data stored in this object.
	 *
	 * Establishes a symmetry mapping between atoms in the system, keeping track of the
	 * set of irreducible atom indices and the information necessary to create the star of atoms.
	 *
	 * @param[in] atoms		A list of the atoms in the system.
	 * @param[in] sym		The symmetry group of the system
	 */
	void initialize(
			std::vector<Atom> const & atoms,
			LatticeStructure::Symmetry const & sym);

	/**
	 * Get the atom index an atom is mapped to by a symmetry operation.
	 *
	 * @param[in] symmetryOpIndex	Index of the symmetry operation that is applied to the atom.
	 * @param[in] atomIndex			Atom index the symmetry operation is applied to.
	 * @return						The index in the original ordering where the atom is mapped to.
	 */
	int atom_rot_map(int symmetryOpIndex, int atomIndex) const;

	/**
	 * Get a list of atom indices that are forming the irreducible set.
	 *
	 * @return	a reference to a list that contains the reducible atom indices that form the irreducible set.
	 */
	std::vector<int> const & get_list_irreducible_atoms() const;

	/**
	 * Get the star of reducible atom indices.
	 *
	 * @param atomIndex		A reducible atom index. Can be any member of the star of indices.
	 * @return				a list of pairs where for any element the first index is a reducible atom index and the second the index of
	 * 						the symmetry operation that rotates the irreducible atom to the reducible one.
	 */
	std::vector<std::pair<int,int>> const & get_star_atom_indices(int atomIndex) const;

	/**
	 * See if an atom is irreducible.
	 * @return	true, if the atom is the reference, irreducible atom, arbitrarily chosen from the star of symmetry related atoms.
	 */
	bool check_atom_is_irreducible(int atomIndex) const;

private:

	/// see symmetry::atom_transform_map()
	std::vector<std::vector<int>> atomSymMap_;

	/// see symmetry::atom_transform_map_inverse()
	std::vector<std::vector<int>> atomSymMapInverse_;

	std::vector<int> irreducibleToReducible_;

	std::vector<Atom> irredAtoms_;

	std::vector<int> redToIrredAtoms_;

	std::vector<int> symRedToIrredAtoms_;

	std::vector< std::vector<int> > irredToRedAtoms_;

	std::vector< std::vector<int> > symIrredToRedAtoms_;

	std::vector<std::vector<std::pair<int,int>>> starAtomIndices_;
};

} /* namespace LatticeStructure */
} /* namespace elephon */

#endif /* ELEPHON_LATTICESTRUCTURE_ATOMSYMMETRYCONNECTION_H_ */
