/*	This file EliashbergSingleRun.cpp is part of elephon.
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
 *  Created on: May 16, 2018
 *      Author: A. Linscheid
 */

#include "EliashbergEquations/EliashbergSingleRun.h"
#include "EliashbergEquations/EliashbergSingleRunData.h"
#include "EliashbergEquations/IsotropicEliashbergEquations.h"
#include "EliashbergEquations/EliashbergModuleConstants.h"

namespace elephon {
namespace EliashbergEquations {

EliashbergSingleRun::EliashbergSingleRun(std::shared_ptr<EliashbergSingleRunData> runData_ptr)
{
	runData_ = runData_ptr;
}

bool
EliashbergSingleRun::solve_Eliashberg()
{
	IsotropicEliashbergEquations eliashbergEquations(runData_);
	bool converged = eliashbergEquations.converge(runData_->get_run_constants()->get_max_number_iterations());
	return converged;
}

EliashbergSingleRun::T
EliashbergSingleRun::get_maximal_gap_magnitude() const
{
	return *std::max_element(runData_->get_gap().get_data_ptr(), runData_->get_gap().get_end_data_ptr());
}

void
EliashbergSingleRun::rescale_gap_function(T newMaxGapValue)
{
	T scale = newMaxGapValue / this->get_maximal_gap_magnitude();
	for (auto it = runData_->access_gap().access_data_ptr(); it != runData_->access_gap().get_end_data_ptr(); ++it)
		*it *= scale;
}

} /* namespace EliashbergEquations */
} /* namespace elephon */
