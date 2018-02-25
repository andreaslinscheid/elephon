/*	This file AtomDisplacementCollection.h is part of elephon.
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
 *  Created on: Feb 10, 2018
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_LATTICESTRUCTURE_ATOMDISPLACEMENTCOLLECTION_H_
#define ELEPHON_LATTICESTRUCTURE_ATOMDISPLACEMENTCOLLECTION_H_

#include "Auxillary/AlignedVector.h"
#include <memory>
#include <vector>

namespace elephon
{
namespace LatticeStructure
{

class UnitCell;
class AtomDisplacement;
class Atom;
class Symmetry;
class LatticeModule;
class AtomSymmetryConnection;

/**
 * Generate and store the set of atomic displacements in the primitive cell.
 */
class AtomDisplacementCollection
{
public:

	/**
	 * Generate displacement patterns in the unit cell.
	 *
	 * @param[in] primitiveCell				The primitive cell
	 * @param[in] displMagn					Magnitude of the displacement
	 * @param[in] symmetricDisplacements	If true, we treat any displacement and the one with opposite direction as equivalent
	 */
	void initialize(
			std::shared_ptr<const UnitCell> primitiveCell,
			bool symmetricDisplacements,
			double displacementMagnitude);

	/**
	 * Build a string with the relative folder structure for a displacement run.
	 *
	 * From a displacement index, resolve the atom index and the irreducible index for this
	 * atom, then return a name appropriate for this setup.
	 *
	 * @param[in] totalIrredIndex	An index in the list of the total number of irreducible displacements.
	 * @return	A string with the name of the irreducible displacement run.
	 */
	std::string get_relative_folder_structure_displ_run(int totalIrredIndex) const;

	/**
	 * Resolve a index in the total list of irreducible displacements into atom and irreducible index at this atom.
	 *
	 * @param totalIrredIndex	An index in the list of the total number of irreducible displacements.
	 * @return	a pair with first being the atom and second being the irreducible displacement index at this atom.
	 */
	std::pair<int,int> get_total_irred_index_to_atom_and_rel_irred_index(int totalIrredIndex) const;

	/**
	 * Get the total number of symmetry inquivalent displacements in the unit cell.
	 * @return the total number of symmetry inquivalent displacements
	 */
	int get_tota_num_irred_displacements() const;

	/**
	 * See if an atom is irreducible.
	 * @return	true, if the atom is the reference, irreducible atom, arbitrarily chosen from the star of symmetry related atoms.
	 */
	bool check_atom_is_irreducible(int atomIndex) const;

	/**
	 * Get the number of irreducible displacements for an atom refered to by the index.
	 *
	 * @param atomIndex	The atom index.
	 * @return	number of irreducible displacements of this atom.
	 */
	int get_num_irred_displacements_for_atom(int atomIndex) const;

	/**
	 * Get number of reducible atoms in the primitive cell.
	 *
	 * @return	the number of reducible atoms in the primitive cell.
	 */
	int get_num_atoms_primitive_cell() const;

	/**
	 * The number of reducible displacements for an atom.
	 * @param atomIndex		Atom index in the primitive cell.
	 * @return				The number of reducible displacements.
	 */
	int get_num_red_displacements_for_atom(int atomIndex) const;

	/**
	 * Get information on how the list of reducible displacement is related to the irreducible ones.
	 *
	 * @param atomIndex			The atom index in the primtive cell
	 * @return					For the ith element in \p reducibleDispl, a pair of ints with the first
	 * 							indicating the index in the list of irreducible displacements and the second
	 * 							the index in the site symmetry group taking one from the irreducible to the reducible
	 * 							displacement.
	 */
	std::vector<std::pair<int,int>> get_symmetry_relation_red_displacements_for_atom(int atomIndex) const;

	/**
	 * Get a list of reducible atomic displacements.
	 *
	 * The relation to the irreducible ones is given with get_symmetry_relation_red_displacements_for_atom()
	 *
	 * @param atomIndex			The atom index in the primtive cell for which we give displacements
	 * @return					The vector of reducible displacements.
	 */
	std::vector<AtomDisplacement> get_red_displacements_for_atom(int atomIndex) const;

	/**
	 * Obtain the irreducible displacements.
	 *
	 * @return A vector of pairs reset to the number of irreducible atoms where the first element of
	 * 		   the pair is the atom index and the second is a vector with displacements belonging to
	 * 		   this atom.
	 */
	std::vector<std::pair<int, std::vector<AtomDisplacement>>> get_irreducible_displacements() const;

	/**
	 * Compute the pseudo inverse of the displacement matrix.
	 *
	 * @param atomIndex		The atom index in the primitive cell.
	 * @param pInv			A multiarray that will be resized to fit the 3 x <number of reducible displacement for this atom>
	 * 						elements containing the pseudo inverse of the displacement matrix.
	 */
	void generate_pseudo_inverse_reducible_displacements_for_atom(
			int atomIndex,
			Auxillary::Multi_array<double,2> & pInv) const;

	/**
	 * Check if displacements of opposite direction are treated equivalent.
	 * @return	return true if displacements of same direction but opposite sign are treated as equivalent.
	 */
	bool get_treat_displacements_pm_symmetric() const;
private:

	std::shared_ptr<const LatticeStructure::AtomSymmetryConnection> atomSym_;

	/// upon generation of displacement, we only store these for symmetry-irreducible atoms
	/// This map tells for a given atom index, the number in the list of displacements.
	std::vector<int> mapReducibleToIrreducibleAtoms_;

	std::vector<std::pair<int,int>> totalIrredToAtomAndRelIrred_;

	std::vector<std::vector<AtomDisplacement>> irredDispl_;

	std::vector<std::vector<AtomDisplacement>> redDispl_;

	std::vector<std::vector<int>> symRedToIrredDispl_;

	std::vector<std::vector<int>> redToIrredDispl_;

	std::vector<std::vector<std::vector<int>>> symIrredToRedDispl_;

	std::vector<std::vector<std::vector<int>>> irredToRedDispl_;

	bool symmetricDisplacements_ = true;

	double displMagn_;

	void generate_displacements(
			std::shared_ptr<const UnitCell> primitiveCell,
			bool symmetricDisplacements,
			double displMagn);

	void get_site_displacements(
			LatticeStructure::LatticeModule const & lattice,
			LatticeStructure::Atom const & atomicSite,
			bool symmetricDisplacements,
			LatticeStructure::Symmetry const & siteSymmetry,
			double displMagn,
			std::vector<AtomDisplacement> & irreducible,
			std::vector<AtomDisplacement> & reducible,
			std::vector<int> & redToIrredDispl,
			std::vector<int> & symRedToIrredDispl,
			std::vector< std::vector<int> > & irredToRedDispl,
			std::vector< std::vector<int> > & symIrredToRedDispl) const;
};

} /* namespace LatticeStructure */
} /* namespace elephon */

#endif /* ELEPHON_LATTICESTRUCTURE_ATOMDISPLACEMENTCOLLECTION_H_ */
