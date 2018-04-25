/*	This file Atom.h is part of elephon.
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

#ifndef ELEPHON_LATTICESTRUCTURE_ATOM_H_
#define ELEPHON_LATTICESTRUCTURE_ATOM_H_

#include "LatticeStructure/Symmetry.h"
#include "symmetry/SymmetryOperation.h"
#include <string>
#include <vector>

namespace elephon
{
namespace LatticeStructure
{

class AtomDisplacement;

/**
 * Basic container of an individual atom.
 */
class Atom
{
public:

	Atom();

	Atom(double mass, std::string kind, std::vector<double> pos,std::vector<bool> frozen, double gridPrec = 1e-6);

	std::string get_kind() const;

	double get_mass() const;

	double get_position_precision() const;

	std::vector<double> const & get_position() const;

	void set_position(std::vector<double> newPos);

	std::vector<bool> get_movement_fixed() const;

	void transform( symmetry::SymmetryOperation const & sop  );

	/**
	 * Induce a translation of the object's position as defined by the AtomDisplacement.
	 *
	 * @param[in] displ	The object definiting the direction and magnitude.
	 */
	void apply_displacement(AtomDisplacement const & displ);

	/**
	 * Establish an ordering of atoms for easier lookup.
	 *
	 * The two atoms compared must share the same grid precision.
	 *
	 * @param a1	one atom.
	 * @param a2	another atom.
	 * @return	true if a1's kind is equal to a2 and they a1's z coordinate is smaller than a2's. If they are equal up to within
	 * 			gridprec, we compare y in the same manner and then x.
	 */
	friend bool operator< (Atom const & a1, Atom const & a2);

private:

	double mass_ = 0.0;

	std::string kind_;

	std::vector<double> pos_;

	std::vector<bool> frozen_;

	double prec_ = 1e-6;

	void map_pos_back_1BZ();
};

/**
 * Compare two atoms if they are equal.
 *
 *	The two atoms compared must share the same grid precision.
 *
 * @param a1	one atom
 * @param a2	another atom
 * @return	true if they are the same kind and share the same position up to within gridprec..
 */
bool operator== (Atom const & a1, Atom const & a2);

} /* namespace LatticeStructure */
} /* namespace elephon */
#endif /* ELEPHON_LATTICESTRUCTURE_ATOM_H_ */
