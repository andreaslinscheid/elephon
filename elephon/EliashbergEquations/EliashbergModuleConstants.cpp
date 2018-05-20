/*	This file EliashbergModuleConstants.cpp is part of elephon.
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

#include "EliashbergEquations/EliashbergModuleConstants.h"

namespace elephon {
namespace EliashbergEquations {

void
EliashbergModuleConstants::initialize(IOMethods::InputOptions const & opts,
		PhononStructure::AlphaSquaredF a2F)
{
	MatsubaraEnergyCutoff_ = opts.get_matsC();
	a2F_ = std::move(a2F);
	isotropicCoulombMuStar_ = opts.get_muStar();
	if (opts.get_muStarC().empty()){
		isotropicCoulombCutoff_ = 10.0*a2F_.get_max_frequency_range();
	} else
	{
		isotropicCoulombCutoff_ = opts.get_muStarC()[0];
	}
	convergenceThreshold_ = opts.get_EliCTh();
	thresholdConsiderZero_ = opts.get_EliThrZero();
	mixingParameter_ = opts.get_EliMix();
	maxNumIterations_ = opts.get_EliNMax();
	shapeConvergenceImprovementThreshold_ = opts.get_EliCImp();
}

PhononStructure::AlphaSquaredF const &
EliashbergModuleConstants::get_isotropic_elphon_coupling() const
{
	return a2F_;
}

double
EliashbergModuleConstants::get_Matsubara_cutoff() const
{
	return MatsubaraEnergyCutoff_;
}

std::vector<double> const &
EliashbergModuleConstants::get_isotropic_coulomb() const
{
	return isotropicCoulombMuStar_;
}

double
EliashbergModuleConstants::get_isotropic_coulomb_cutoff() const
{
	return isotropicCoulombCutoff_;
}

int
EliashbergModuleConstants::get_number_bands() const
{
	return a2F_.get_num_bands();
}

double
EliashbergModuleConstants::get_convergence_threshold() const
{
	return convergenceThreshold_;
}

double
EliashbergModuleConstants::get_considered_zero_threshold() const
{
	return thresholdConsiderZero_;
}

int
EliashbergModuleConstants::get_max_number_iterations() const
{
	return maxNumIterations_;
}

double
EliashbergModuleConstants::get_mixing_parameter() const
{
	return mixingParameter_;
}

double
EliashbergModuleConstants::get_shapeConvergence_improvement_threshold() const
{
	return shapeConvergenceImprovementThreshold_;
}

} /* namespace EliashbergEquations */
} /* namespace elephon */
