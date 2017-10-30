/*	This file LocalDensityOfStates.h is part of elephon.
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
 *  Created on: Jul 24, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_ELECTRONICSTRUCTURE_LOCALDENSITYOFSTATES_H_
#define ELEPHON_ELECTRONICSTRUCTURE_LOCALDENSITYOFSTATES_H_

#include "ElectronicStructure/Wavefunctions.h"
#include "ElectronicStructure/GradientFFTReciprocalGrid.h"
#include "IOMethods/ResourceHandler.h"
#include "ElectronicStructure/TetrahedraIsosurface.h"
#include <vector>
#include <string>
#include <memory>

namespace elephon
{
namespace ElectronicStructure
{

class LocalDensityOfStates
{
public:

	void compute_ldos(
			std::vector<double> energies,
			std::shared_ptr<const Wavefunctions> wfcts,
			std::shared_ptr<const LatticeStructure::UnitCell> unitcell,
			int nkpointsPerSurface,
			std::shared_ptr<const LatticeStructure::RegularBareGrid> realSpaceRes,
			std::shared_ptr<const TetrahedraIsosurface> tetraIso,
			bool symmetrize);

	/**
	 * Compute the local density of states.
	 *
	 * This version creates local variables with wave functions, fermi surface, grids, ect.
	 *
	 * @param energies	Compute the the LDOS for each of these energies relative
	 * 					to the Fermi level.
	 * @param loader	The electronic structure code interface through which we
	 * 					load all needed data.
	 */
	void compute_ldos(
			std::vector<double> const & energies,
			std::shared_ptr<IOMethods::ResourceHandler> loader );

	void write_file( std::string const & filename, bool binary = false) const;
private:

	std::shared_ptr<const LatticeStructure::RegularBareGrid> rsDims_;

	std::vector<double> isoEnergies_;

	std::shared_ptr<const LatticeStructure::UnitCell> uc_;

	std::vector<double> ldos_;
};

} /* namespace ElectronicStructure */
} /* namespace elephon */

#endif /* ELEPHON_ELECTRONICSTRUCTURE_LOCALDENSITYOFSTATES_H_ */
