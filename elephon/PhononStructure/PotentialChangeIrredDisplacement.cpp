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
#include "AtomicSite/SphericalHarmonicExpansion.h"
#include "LatticeStructure/AtomDisplacement.h"
#include "Algorithms/rotation_matrix_connecting_unit_vectors.h"

namespace elephon
{
namespace PhononStructure
{

void
PotentialChangeIrredDisplacement::initialize(
		std::shared_ptr<const LatticeStructure::AtomDisplacement> displ,
		std::vector<double> const & regularGridGroundStatePotential,
		std::vector<std::shared_ptr<const AtomicSite::SphericalHarmonicExpansion>> radialGroundStatePotential,
		std::vector<double> const & regularGridDisplacedPotential,
		std::vector<std::shared_ptr<const AtomicSite::SphericalHarmonicExpansion>> radialDisplacedPotential,
		std::shared_ptr<const LatticeStructure::RegularBareGrid> unitcellGrid,
		std::shared_ptr<const LatticeStructure::RegularBareGrid> supercellGrid )
{
	assert(regularGridGroundStatePotential.size() == unitcellGrid->get_num_points());
	assert(regularGridDisplacedPotential.size() == supercellGrid->get_num_points());
	// check that both ways of specifying the number of primitive cells in a supercell is equal
	assert(radialDisplacedPotential.size()/radialGroundStatePotential.size()
			== supercellGrid->get_num_points()/unitcellGrid->get_num_points());

	// regular grid data is first
	auto ucDim = unitcellGrid->get_grid_dim();
	std::vector<int> xyz(3);
	regularGridData_.resize(regularGridDisplacedPotential.size());
	for (int irsc = 0 ; irsc < regularGridDisplacedPotential.size(); ++irsc)
	{
		supercellGrid->get_reducible_to_xyz(irsc, xyz);
		for (int i = 0 ; i < 3 ; ++i)
			xyz[i] = xyz[i]%ucDim[i];
		auto irUC = unitcellGrid->get_xyz_to_reducible(xyz);
		assert((irUC >=0) && (irUC<regularGridGroundStatePotential.size()));
		regularGridData_[irsc] = regularGridGroundStatePotential[irUC]-regularGridDisplacedPotential[irsc];
	}

	// identify the connection of the primitive unit cell to the supercell where one atom is displaced

	// with this connection map of atomic sites, evaluate the potential difference on the position of the unperturbed grid
	// For all but the displaced radial position, this is straight forward. We exempt this site from the calculation until further.
	radialGridData_.reserve(radialDisplacedPotential.size());

	// perform the calculation of the displaced site. The formula is given in the documentation of the main class.
}

} /* namespace PhononStructure */
} /* namespace elephon */
