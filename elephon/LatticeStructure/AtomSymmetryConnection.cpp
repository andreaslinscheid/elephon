/*	This file AtomSymmetryConnection.cpp is part of elephon.
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

#include "LatticeStructure/AtomSymmetryConnection.h"
#include "symmetry/atom_transform_map.h"
#include "LatticeStructure/SymmetryReduction.h"

namespace elephon
{
namespace LatticeStructure
{

void
AtomSymmetryConnection::initialize(
		std::vector<Atom> const & atoms,
		LatticeStructure::Symmetry const & sym)
{
	symmetry::atom_transform_map(atoms, sym, atomSymMap_);
	symmetry::atom_transform_map_inverse(atoms, sym, atomSymMapInverse_);

	SymmetryReduction<Atom>(
			sym,
			atoms,  irredAtoms_,
			redToIrredAtoms_, symRedToIrredAtoms_,
			irredToRedAtoms_, symIrredToRedAtoms_);

	irreducibleToReducible_.resize(irredToRedAtoms_.size(), -1);
	for (int iaIrred = 0 ; iaIrred < irredToRedAtoms_.size(); ++iaIrred)
	{
		for (int iStar = 0 ; iStar < irredToRedAtoms_[iaIrred].size(); ++iStar)
		{
			if( symIrredToRedAtoms_[iaIrred][iStar] == sym.get_identity_index() )
			{
				irreducibleToReducible_[iaIrred] = irredToRedAtoms_[iaIrred][iStar];
				break;
			}
		}
	}
	assert(*std::min_element(irreducibleToReducible_.begin(), irreducibleToReducible_.end()) >= 0);

	starAtomIndices_.resize(irredToRedAtoms_.size());
	for (int iaIrred = 0 ; iaIrred < irredToRedAtoms_.size(); ++iaIrred)
	{
		starAtomIndices_[iaIrred].resize(irredToRedAtoms_[iaIrred].size());
		for (int iStar = 0 ; iStar < irredToRedAtoms_[iaIrred].size(); ++iStar)
		{
			starAtomIndices_[iaIrred][iStar] = std::make_pair(symIrredToRedAtoms_[iaIrred][iStar], irredToRedAtoms_[iaIrred][iStar]);
		}
	}
}

int
AtomSymmetryConnection::atom_rot_map(int symmetryOpIndex, int atomIndex) const
{
	assert((atomSymMapInverse_.size()>symmetryOpIndex) && (symmetryOpIndex>=0));
	assert((atomSymMapInverse_[symmetryOpIndex].size()>atomIndex) && (atomIndex>=0));
	return atomSymMapInverse_[symmetryOpIndex][atomIndex];
}

std::vector<int> const &
AtomSymmetryConnection::get_list_irreducible_atoms() const
{
	return irreducibleToReducible_;
}

std::vector<std::pair<int,int>> const &
AtomSymmetryConnection::get_star_atom_indices(int atomIndex) const
{
	assert((atomIndex>=0)&&(atomIndex<redToIrredAtoms_.size()));
	int irreducibleIndex = redToIrredAtoms_[atomIndex];
	return starAtomIndices_[irreducibleIndex];
}

bool
AtomSymmetryConnection::check_atom_is_irreducible(int atomIndex) const
{
	assert((atomIndex>=0)&&(atomIndex<redToIrredAtoms_.size()));
	int irreducibleIndex = redToIrredAtoms_[atomIndex];
	assert((irreducibleIndex>=0)&&(irreducibleIndex<irreducibleToReducible_.size()));
	return irreducibleToReducible_[irreducibleIndex] == atomIndex;
}

} /* namespace LatticeStructure */
} /* namespace elephon */
