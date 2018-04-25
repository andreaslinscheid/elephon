/*	This file Forces.h is part of elephon.
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
 *  Created on: Feb 23, 2018
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_PHONONSTRUCTURE_FORCES_H_
#define ELEPHON_PHONONSTRUCTURE_FORCES_H_

#include "Auxillary/AlignedVector.h"
#include <memory>

namespace elephon
{
namespace LatticeStructure { class AtomDisplacementCollection; };
namespace IOMethods { class ElectronicStructureCodeInterface; };
namespace LatticeStructure { class Symmetry; };
namespace LatticeStructure { class Atom; };

namespace PhononStructure
{

/**
 * A class that manages the forces as read in from the electronic structure calculator.
 */
class Forces
{
public:

	void initialize(
			std::shared_ptr<const LatticeStructure::AtomDisplacementCollection> displ,
			std::shared_ptr<IOMethods::ElectronicStructureCodeInterface> loader );

	/**
	 * Obtain the irreducible set of forces induces by a displacement of an atom in the primitve cell.
	 *
	 * @param[in] atomIndex			The index of the atom in the primitive cell.
	 * @param[in] irredDisplIndex	The index of the irreducible displacement at this atom.
	 * @return					A reference to a \f$ N_{{\rm Atoms-SC}}\times 3\f$ multiarray for each atom
	 * 							in the supercell the 3 component force in cartesian coordinates.
	 */
	Auxillary::Multi_array<double,2> const & get_forces_for_atom( int atomIndex, int irredDisplIndex) const;

	int get_num_total_irred_displacements() const;

	/**
	 *	Obtain the forces induced by all reducible displacements.
	 *
	 *	The data is generated from the one for the irreducible dispalements by applying the symmetry operations.
	 *
	 * @param[in] siteSymmetry			The site symmetry group for the atom with the index \p primitiveAtomIndex
	 * @param[in] primitiveAtomIndex	Index of the primitive atom that is displaced.
	 * @param[in] supercellAtomIndex	Index of the same atom as \p primitiveAtomIndex but in the listing of the supercell.
	 * @param[in] supercellAtomsList	A list of atoms in the supercell.
	 * @param[out] symExpandedForces	Resized to fit NReducible x NatomsSC x 3, i.e. for each reducible displacement
	 * 									in the same order as provided by the underlying LatticeStructure::AtomDisplacementCollection,
	 * 									the forces induced on the atoms in the supercell (by their index) in cartesian
	 * 									x,y and z as the last direction
	 */
	void site_symmetry_expand_data(
			LatticeStructure::Symmetry const & siteSymmetry,
			int primitiveAtomIndex,
			int supercellAtomIndex,
			std::vector<LatticeStructure::Atom> supercellAtomsList,
			Auxillary::Multi_array<double,3> & symExpandedForces ) const;

private:

	int totalNumIrredDispl_;

	Auxillary::Multi_array<int,2> indexMap_;

	std::vector<Auxillary::Multi_array<double,2>> forceData_;

	std::shared_ptr<const LatticeStructure::AtomDisplacementCollection> displColl_;
};

} /* namespace PhononStructure */
} /* namespace elephon */

#endif /* ELEPHON_PHONONSTRUCTURE_FORCES_H_ */
