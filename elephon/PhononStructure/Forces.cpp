/*	This file Forces.cpp is part of elephon.
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

#include "PhononStructure/Forces.h"
#include "LatticeStructure/AtomDisplacementCollection.h"
#include "IOMethods/ElectronicStructureCodeInterface.h"
#include "boost/filesystem.hpp"

namespace elephon
{
namespace PhononStructure
{

void
Forces::initialize(
		std::shared_ptr<const LatticeStructure::AtomDisplacementCollection> displ,
		std::shared_ptr<IOMethods::ElectronicStructureCodeInterface> loader )
{
	auto phononDir = boost::filesystem::path( loader->get_optns().get_elphd());

	// set up the map
	int maxNumIrred = 0;
	const int numAtomPC = displ->get_num_atoms_primitive_cell();
	for (int iAPC = 0 ; iAPC < numAtomPC ; ++iAPC)
		maxNumIrred = std::max(maxNumIrred, displ->get_num_irred_displacements_for_atom(iAPC));

	indexMap_.resize(boost::extents[numAtomPC][maxNumIrred]);
	std::fill_n(indexMap_.data(), indexMap_.size(), -1);

	// read the data
	totalNumIrredDispl_ = displ->get_tota_num_irred_displacements();
	forceData_.resize(totalNumIrredDispl_);
	std::vector<double> thisForces;
	int numAtomSC = -1;
	for ( int idispl = 0 ; idispl < totalNumIrredDispl_; ++idispl )
	{
		std::string dir = (phononDir / displ->get_relative_folder_structure_displ_run(idispl)).string();
		loader->read_forces(dir, thisForces);
		if ( (numAtomSC >= 0) and (numAtomSC != thisForces.size()/3))
			throw std::runtime_error(std::string()+
					"Problem found that the number of atoms deduced from the forces is not equal for all displacements\n"
					"Problem occured for displacement dir: "+dir);
		numAtomSC = thisForces.size()/3;

		forceData_[idispl].resize(boost::extents[numAtomSC][3]);
		std::copy(thisForces.begin(), thisForces.end(), forceData_[idispl].data());

		auto atomAndIrredIndex = displ->get_total_irred_index_to_atom_and_rel_irred_index(idispl);
		assert(atomAndIrredIndex.first<indexMap_.shape()[0]);
		indexMap_[atomAndIrredIndex.first][atomAndIrredIndex.second] = idispl;
	}

	displColl_ = displ;
}

Auxillary::Multi_array<double,2> const &
Forces::get_forces_for_atom( int atomIndex, int irredDisplIndex) const
{
	assert((indexMap_[atomIndex][irredDisplIndex] >= 0) && (indexMap_[atomIndex][irredDisplIndex] < forceData_.size()));
	return forceData_[indexMap_[atomIndex][irredDisplIndex]];
}

void
Forces::site_symmetry_expand_data(LatticeStructure::Symmetry const & siteSymmetry,
		 int primitiveAtomIndex, int supercellAtomIndex,
		std::vector<LatticeStructure::Atom> supercellAtomsList,
		Auxillary::Multi_array<double,3> & symExpandedForces ) const
{
	const int nRedDispl = displColl_->get_num_red_displacements_for_atom(primitiveAtomIndex);
	const int nASC = supercellAtomsList.size();
	symExpandedForces.resize(boost::extents[nRedDispl][nASC][3]);

	// obtain a symmetry mapping of the local environment at atom 'primitiveAtomIndex'
	// a.k.a 'supercellAtomIndex' that tells for each given symmetry operation where
	// an atom in the supercell is mapped to.
	assert((supercellAtomIndex>=0) && (supercellAtomIndex<supercellAtomsList.size()));
	auto const centerAtomPosition = supercellAtomsList[supercellAtomIndex].get_position();
	std::map<LatticeStructure::Atom,int> shiftedSCLookup;
	for (int iASC = 0 ; iASC < nASC; ++iASC)
	{
		auto pos = supercellAtomsList[iASC].get_position();
		for (int xi = 0 ; xi < 3; ++xi)
			pos[xi] -= centerAtomPosition[xi];
		supercellAtomsList[iASC].set_position(std::move(pos));
		shiftedSCLookup.insert(std::make_pair(supercellAtomsList[iASC],iASC));
	}

	auto starDisplacements = displColl_->get_symmetry_relation_red_displacements_for_atom(primitiveAtomIndex);
	for (int redDispl = 0 ; redDispl < starDisplacements.size() ; ++redDispl  )
	{
		const int irredDispl = starDisplacements[redDispl].first;
		const int symIrredToRedDispl = starDisplacements[redDispl].second;

		// in a first step, we copy the force to the right atom location, as taken to by the symmetry
		// and in the second step, we transform the entire set in x,y and z using the symmetry operation
		auto sop = siteSymmetry.get_sym_op(symIrredToRedDispl);
		Auxillary::Multi_array<double,2> const & irredForce = this->get_forces_for_atom(primitiveAtomIndex,irredDispl);

		for (int iASC = 0 ; iASC < nASC; ++iASC)
		{
			auto a_rot = supercellAtomsList[iASC];
			a_rot.transform(sop);
			auto it = shiftedSCLookup.find(a_rot);
			if ( it == shiftedSCLookup.end() )
				throw std::logic_error(std::string("Problem locating supercell atom ")+std::to_string(iASC)+
						" from the point of view of primitive cell atom "+std::to_string(primitiveAtomIndex)+
						" for sym-op no. " + std::to_string(symIrredToRedDispl));
			assert((it->second>=0) && (it->second<nASC));
			for (int xi = 0 ; xi < 3 ; ++xi)
				symExpandedForces[redDispl][it->second][xi] = irredForce[iASC][xi];
		}

		siteSymmetry.rotate_cartesian(symIrredToRedDispl, &symExpandedForces[redDispl][0][0], &symExpandedForces[redDispl][0][0]+3*nASC);
	}
}

int
Forces::get_num_total_irred_displacements() const
{
	return totalNumIrredDispl_;
}

} /* namespace PhononStructure */
} /* namespace elephon */
