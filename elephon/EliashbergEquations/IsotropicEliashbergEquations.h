/*	This file IsotropicEliashbergEquations.h is part of elephon.
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

#ifndef ELEPHON_ELIASHBERGEQUATIONS_ISOTROPICELIASHBERGEQUATIONS_H_
#define ELEPHON_ELIASHBERGEQUATIONS_ISOTROPICELIASHBERGEQUATIONS_H_

#include <memory>

namespace elephon {
namespace EliashbergEquations {

class EliashbergSingleRunData;
class EliashbergMixingModule;

class IsotropicEliashbergEquations
{
public:
	IsotropicEliashbergEquations(std::shared_ptr<EliashbergSingleRunData> runData);

	bool converge(int maximalNumberOfIterations);
private:
	std::shared_ptr<EliashbergSingleRunData> runData_;

	std::shared_ptr<EliashbergMixingModule> mixing_;

	bool check_convergence() const;

	void gap_analysis_and_convergence_improvement();

	void mix_iterations();

	void gap_equation();

	void frequRenorm_equation();
};

} /* namespace EliashbergEquations */
} /* namespace elephon */

#endif /* ELEPHON_ELIASHBERGEQUATIONS_ISOTROPICELIASHBERGEQUATIONS_H_ */
