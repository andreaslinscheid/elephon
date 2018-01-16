/*	This file atom_transform_map.cpp is part of elephon.
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
 *  Created on: Jan 15, 2018
 *      Author: A. Linscheid
 */

#include "symmetry/atom_transform_map.h"
#include <map>

namespace elephon
{
namespace symmetry
{

void
atom_transform_map_both(
		std::vector<LatticeStructure::Atom> const & atoms,
		LatticeStructure::Symmetry const & symmetry,
		std::vector< std::vector<int> > & rotAtomsMap,
		bool inverse)
{
	int iA = atoms.size();
	int iS = symmetry.get_num_symmetries();
	//We create a lookup for atoms, then we apply the rotation and
	//try to discover the transformed position in the set
	std::map<LatticeStructure::Atom, int> loopup;
	for ( int i = 0 ; i < iA ; ++i)
		loopup.insert( std::move(std::make_pair(atoms[i],i)) );

	rotAtomsMap = std::vector< std::vector<int> >(iS, std::vector<int>(iA) );
	for ( int isym = 0 ; isym < iS; ++isym)
	{
		//rotate all atoms
		auto rotAtoms = atoms;
		for ( auto &a : rotAtoms )
			a.transform(symmetry.get_sym_op(isym));

		for ( int i = 0 ; i < iA ; ++i )
		{
			auto it = loopup.find( rotAtoms[i] );
			if ( it == loopup.end() )
				throw std::logic_error("The set of atoms is not closed "
						"under symmetry operations which can't be.");
			if ( inverse )
				rotAtomsMap[isym][it->second] = i;
			else
				rotAtomsMap[isym][i] = it->second;
		}
	}
}

void
atom_transform_map(
		std::vector<LatticeStructure::Atom> const & atoms,
		LatticeStructure::Symmetry const & symmetry,
		std::vector< std::vector<int> > & rotAtomsMap)
{
	atom_transform_map_both(atoms, symmetry, rotAtomsMap, false);
}

void
atom_transform_map_inverse(
		std::vector<LatticeStructure::Atom> const & atoms,
		LatticeStructure::Symmetry const & symmetry,
		std::vector< std::vector<int> > & rotAtomsMap)
{
	atom_transform_map_both(atoms, symmetry, rotAtomsMap, true);
}

} /* namespace Symmetry */
} /* namespace elephon */
