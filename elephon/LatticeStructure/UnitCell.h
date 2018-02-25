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
#include <memory>

namespace elephon
{
namespace LatticeStructure
{

class AtomSymmetryConnection;

/**
 * A collection of data and methods that describes the physical system in a periodic lattice.
 *
 */
class UnitCell
{
public:
	/**
	 * Define the system.
	 *
	 * @param atoms		List of atoms in the unit cell specifing the basis of the lattice.
	 * @param lattice	A lattice structure specifying the cell.
	 * @param sym		The symmetry of the system.
	 */
	void initialize(
			std::vector<LatticeStructure::Atom> atoms,
			LatticeStructure::LatticeModule lattice,
			LatticeStructure::Symmetry sym);

	/**
	 * Construct a supercell from the primitive unit cell.
	 *
	 * @todo refactor to independent class that also allows non-trivial supercells.
	 *
	 * @param[in] scx	The multiple of the x dimension
	 * @param[in] scy	The multiple of the y dimension
	 * @param[in] scz	The multiple of the z dimension
	 * @return	An initialized super-unitcell.
	 */
	UnitCell build_supercell(int scx, int scy, int scz) const;

	void displace_atom( AtomDisplacement const& displ );

	void get_site_displacements(Atom const & atomicSite,
			bool symmetricDisplacements,
			LatticeStructure::Symmetry const & siteSymmetry,
			double displMagn,
			std::vector<AtomDisplacement> & irreducible,
			std::vector<AtomDisplacement> & reducible,
			std::vector<int> & redToIrredDispl,
			std::vector<int> & symRedToIrredDispl,
			std::vector< std::vector<int> > & irredToRedDispl,
			std::vector< std::vector<int> > & symIrredToRedDispl) const;

	std::vector<LatticeStructure::Atom> const & get_atoms_list() const;

	LatticeStructure::LatticeModule const & get_lattice() const;

	double get_alat() const;

	LatticeStructure::Symmetry const & get_symmetry() const;

	LatticeStructure::Symmetry const & get_site_symmetry(int atomIndex) const;

	void set_symmetry_to_lattice(LatticeStructure::Symmetry & symmetry) const;

	void compute_supercell_dim(std::shared_ptr<const UnitCell> supercell, std::vector<int> & supercellDim ) const;

	std::shared_ptr<const LatticeStructure::AtomSymmetryConnection> get_atom_symmetry() const;
private:

	LatticeStructure::LatticeModule lattice_;

	LatticeStructure::Symmetry symmetry_;

	std::vector<LatticeStructure::Atom> atoms_;

	std::vector<LatticeStructure::Symmetry> siteSymmetries_;

	std::shared_ptr<LatticeStructure::AtomSymmetryConnection> atomSymmetryModule_;

	void add_displacement( std::vector<double> direction,
			std::vector<double> const & position,
			bool symmetricDispl,
			std::string const & atomName,
			double gridPrec,
			double magnInAngstroem,
			std::vector<AtomDisplacement> & addtothis) const;

	void generate_site_symmetries(std::vector<LatticeStructure::Atom> const & atoms,
			LatticeStructure::Symmetry const & symmetry,
			std::vector<LatticeStructure::Symmetry> & siteSymmetries) const;
};

} /* namespace LatticeStructure */
} /* namespace elephon */

#endif /* ELEPHON_LATTICESTRUCTURE_UNITCELL_H_ */
