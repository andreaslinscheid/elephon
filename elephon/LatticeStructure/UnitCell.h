/*	This file UnitCell.h is part of elephon.
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
 *  Created on: May 15, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_LATTICESTRUCTURE_UNITCELL_H_
#define ELEPHON_LATTICESTRUCTURE_UNITCELL_H_

#include "LatticeStructure/LatticeModule.h"
#include "LatticeStructure/AtomDisplacement.h"
#include "LatticeStructure/Symmetry.h"
#include "LatticeStructure/Atom.h"

namespace elephon
{
namespace LatticeStructure
{

class UnitCell
{
public:
	void initialize(
			std::vector<LatticeStructure::Atom> atoms,
			LatticeStructure::LatticeModule lattice,
			LatticeStructure::Symmetry sym);

	UnitCell build_supercell(int scx, int scy, int scz) const;

	void displace_atom( AtomDisplacement const& displ );

	void generate_displacements( double displMagn,
			bool symmetricDisplacements,
			std::vector<AtomDisplacement> & reducibleDisplacements,
			std::vector<AtomDisplacement> & irreducibleDisplacements,
			std::vector<int> & redToIrred,
			std::vector<int> & symRedToIrred,
			std::vector< std::vector<int> > & irredToRed,
			std::vector< std::vector<int> > & symIrredToRed) const;

	std::vector<LatticeStructure::Atom> const & get_atoms_list() const;

	std::vector<double> get_lattice_matrix() const;

	double get_alat() const;
private:

	LatticeStructure::LatticeModule lattice_;

	LatticeStructure::Symmetry symmetry_;

	std::vector<LatticeStructure::Atom> atoms_;

	void reduce_symmetry_to_lattice();
};

} /* namespace LatticeStructure */
} /* namespace elephon */

#endif /* ELEPHON_LATTICESTRUCTURE_UNITCELL_H_ */
