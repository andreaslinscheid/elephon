/*	This file EliashbergModuleConstants.h is part of elephon.
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
 *  Created on: May 15, 2018
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_ELIASHBERGEQUATIONS_ELIASHBERGMODULECONSTANTS_H_
#define ELEPHON_ELIASHBERGEQUATIONS_ELIASHBERGMODULECONSTANTS_H_

#include "EliashbergEquations/IsotropicElectronPhononCoupling.h"
#include "PhononStructure/AlphaSquaredF.h"
#include "IOMethods/InputOptions.h"

namespace elephon
{
namespace EliashbergEquations
{
/**
 * A class that stores data which remains constant through out the work done in the Eliashberg Module.
 */
class EliashbergModuleConstants
{
public:

	void initialize( IOMethods::InputOptions const & opts,
			PhononStructure::AlphaSquaredF a2F);

	double get_Matsubara_cutoff() const;

	std::vector<double> const & get_isotropic_coulomb() const;

	double get_isotropic_coulomb_cutoff() const;

	int get_number_bands() const;

	double get_convergence_threshold() const;

	double get_considered_zero_threshold() const;

	int get_max_number_iterations() const;

	double get_mixing_parameter() const;

	PhononStructure::AlphaSquaredF const & get_isotropic_elphon_coupling() const;

	double get_shapeConvergence_improvement_threshold() const;
private:

	double MatsubaraEnergyCutoff_;

	std::vector<double> isotropicCoulombMuStar_;

	double isotropicCoulombCutoff_;

	double convergenceThreshold_;

	double thresholdConsiderZero_;

	double mixingParameter_;

	int maxNumIterations_;

	double shapeConvergenceImprovementThreshold_;

	PhononStructure::AlphaSquaredF a2F_;
};

} /* namespace EliashbergEquations */
} /* namespace elephon */

#endif /* ELEPHON_ELIASHBERGEQUATIONS_ELIASHBERGMODULECONSTANTS_H_ */
