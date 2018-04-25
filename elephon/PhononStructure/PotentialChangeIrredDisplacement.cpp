/*	This file PotentialChangeIrredDisplacement.cpp is part of elephon.
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
 *  Created on: Jan 17, 2018
 *      Author: A. Linscheid
 */

#include "PhononStructure/PotentialChangeIrredDisplacement.h"
#include "LatticeStructure/RegularBareGrid.h"
#include "AtomicSite/AtomSiteData.h"
#include "LatticeStructure/AtomDisplacement.h"
#include "LatticeStructure/Atom.h"
#include "LatticeStructure/PrimitiveToSupercellConnection.h"
#include "Algorithms/rotation_matrix_connecting_unit_vectors.h"
#include "AtomicSite/AtomSiteData.h"
#include <limits>
#include <boost/multi_array.hpp>

namespace elephon
{
namespace PhononStructure
{

void
PotentialChangeIrredDisplacement::initialize(
		LatticeStructure::AtomDisplacement const & atomDispl,
		LatticeStructure::DataRegularAndRadialGrid<double> const & groundStatePotential,
		LatticeStructure::DataRegularAndRadialGrid<double> const & displacedPotential,
		std::shared_ptr<const LatticeStructure::RegularBareGrid> unitcellGrid,
		std::shared_ptr<const LatticeStructure::RegularBareGrid> supercellGrid,
		std::shared_ptr<const LatticeStructure::PrimitiveToSupercellConnection> scPrimMap)
{
#ifndef NDEBUG
	const int nGridRegularGS = std::distance(groundStatePotential.begin_regular_data(), groundStatePotential.end_regular_data());
	const int nGridRegularDS = std::distance(displacedPotential.begin_regular_data(), displacedPotential.end_regular_data());
	assert(nGridRegularGS == unitcellGrid->get_num_points());
	assert(nGridRegularDS == supercellGrid->get_num_points());
	// check that both ways of specifying the number of primitive cells in a supercell is equal
	assert(nGridRegularDS/nGridRegularGS
			== supercellGrid->get_num_points()/unitcellGrid->get_num_points());
#endif

	// regular grid data is first
	auto ucDim = unitcellGrid->get_grid_dim();
	std::vector<int> xyz(3);
	Auxillary::alignedvector::DV regularGridData(supercellGrid->get_num_points());
	for (int irsc = 0 ; irsc < supercellGrid->get_num_points(); ++irsc)
	{
		supercellGrid->get_reducible_to_xyz(irsc, xyz);
		for (int i = 0 ; i < 3 ; ++i)
			xyz[i] = xyz[i]%ucDim[i];
		auto irUC = unitcellGrid->get_xyz_to_reducible(xyz);
		assert((irUC >=0) && (irUC<unitcellGrid->get_num_points()));
		regularGridData[irsc] = *(groundStatePotential.begin_regular_data() +irUC)
								- *(displacedPotential.begin_regular_data() + irsc);
	}

	// build the unperturbed data in the supercell
	std::vector<AtomicSite::AtomSiteData> radialGridData(displacedPotential.get_num_atom_data_sets(), AtomicSite::AtomSiteData());
	Auxillary::Multi_array<int,2> Rvectors;
	Auxillary::Multi_array<int,3> indexMapRvectors;
	scPrimMap->get_supercell_vectors(Rvectors, indexMapRvectors);

	for (int iAPC = 0; iAPC < groundStatePotential.get_num_atom_data_sets(); ++iAPC)
	{
		AtomicSite::AtomSiteData const & d = groundStatePotential.view_radial_data_set(iAPC);
		auto atom = d.get_atom();
		for (int iNumPrim = 0; iNumPrim < Rvectors.shape()[0]; ++iNumPrim)
		{
			AtomicSite::AtomSiteData newData;
			// get the position, add the lattice vector and scale from supercell to primitive cell coordinates
			assert(Rvectors.shape()[1] == 3);
			auto pos = atom.get_position();
			for (int i = 0; i < 3; ++i)
				pos[i] += Rvectors[iNumPrim][i];
			scPrimMap->primitive_to_supercell_coordinates(pos);
			auto atomCopy = atom;
			atomCopy.set_position(pos);
			int atomIndexSupercell = scPrimMap->find_atom(atomCopy, false);
			assert(atomIndexSupercell>=0);
			newData.initialize(std::move(atomCopy), d.get_data());
			radialGridData[atomIndexSupercell] = std::move(newData);
		}
	}

	// reconstruct the displaced atom from the displacement and find it in the unit cell
	LatticeStructure::Atom dACheck(1.0, atomDispl.get_kind(), atomDispl.get_position(), {false, false, false}, atomDispl.get_prec() );
	int displacedAtomIndexPrimitive = scPrimMap->find_atom(dACheck, true);
	if(displacedAtomIndexPrimitive < 0)
		throw std::logic_error("Atom referenced by displ not found in the primitive cell");
	int displacedAtomIndexSupercell = scPrimMap->primitive_to_supercell_atom_index(displacedAtomIndexPrimitive);

	// this calculates the finite displacement coefficient wise for all but the displaced atom
	const AtomicSite::AtomSiteData * dispAtomData_ptr = nullptr;
	bool oneMatched = false;
	for (int iASC = 0; iASC < displacedPotential.get_num_atom_data_sets(); ++iASC)
	{
		AtomicSite::AtomSiteData const & perturbedData = displacedPotential.view_radial_data_set(iASC);
		int atomIndexSupercell = scPrimMap->find_atom(perturbedData.get_atom(), false);
		if ( (atomIndexSupercell < 0) or (atomIndexSupercell == displacedAtomIndexSupercell) )
		{
			dispAtomData_ptr = &perturbedData;
			if( oneMatched )
				throw std::logic_error("Atoms in radialDisplacedPotential overlap or only one atom must be displaced.");
			oneMatched = true;
			continue;
		}
		auto itUnperturb = radialGridData[atomIndexSupercell].edit_data().begin();
		auto itPerturb = perturbedData.get_data().begin();
		for ( ; itUnperturb != radialGridData[atomIndexSupercell].edit_data().end(); ++itUnperturb, ++itPerturb )
		{
			*itUnperturb -= (*itPerturb);
		}
	}
	assert(dispAtomData_ptr != nullptr);

	// perform the calculation of the displaced site. The formula is given in the documentation of the main class.
	// copy the grid but then set it back to the original position - this is where we will evaluate the data
	auto rGridUperturbedPos = dispAtomData_ptr->get_data().get_radial_grid();
	rGridUperturbedPos.set_center(radialGridData[displacedAtomIndexSupercell].get_atom().get_position());

	AtomicSite::SphericalHarmonicExpansion shiftedBackFit;
	shiftedBackFit.fit_to_data(
			dispAtomData_ptr->get_data(),
			dispAtomData_ptr->get_data().get_l_max(),
			std::move(rGridUperturbedPos));
	AtomicSite::AtomSiteData newDataDisplacedSite;
	newDataDisplacedSite.initialize(dispAtomData_ptr->get_atom(), std::move(shiftedBackFit));

	// take the difference also for the displaced atom and the -now- shifted data.
	auto itUnperturb = radialGridData[displacedAtomIndexSupercell].edit_data().begin();
	auto itPerturb = newDataDisplacedSite.get_data().begin();
	for ( ; itUnperturb != radialGridData[displacedAtomIndexSupercell].edit_data().end(); ++itUnperturb, ++itPerturb )
	{
		*itUnperturb -= (*itPerturb);
	}

	data_.initialize(*supercellGrid, std::move(regularGridData), std::move(radialGridData));
}

int
PotentialChangeIrredDisplacement::get_max_num_radial_elements() const
{
	return data_.get_max_num_radial_elements();
}

int
PotentialChangeIrredDisplacement::get_max_angular_moment() const
{
	return data_.get_max_angular_moment();
}

int
PotentialChangeIrredDisplacement::get_max_num_angular_moment_channels() const
{
	return data_.get_max_num_angular_moment_channels();
}

void
PotentialChangeIrredDisplacement::transform(symmetry::SymmetryOperation const & sop)
{
	data_.transform(sop);
}

Auxillary::alignedvector::DV::const_iterator
PotentialChangeIrredDisplacement::begin_regular_data() const
{
	return data_.begin_regular_data();
}

Auxillary::alignedvector::DV::const_iterator
PotentialChangeIrredDisplacement::end_regular_data() const
{
	return data_.end_regular_data();
}

Auxillary::alignedvector::ZV::const_iterator
PotentialChangeIrredDisplacement::begin_radial_data(int atomIndex) const
{
	return data_.begin_radial_data(atomIndex);
}

Auxillary::alignedvector::ZV::const_iterator
PotentialChangeIrredDisplacement::end_radial_data(int atomIndex) const
{
	return data_.end_radial_data(atomIndex);
}

LatticeStructure::DataRegularAndRadialGrid<double> const &
PotentialChangeIrredDisplacement::view_data() const
{
	return data_;
}

} /* namespace PhononStructure */
} /* namespace elephon */
