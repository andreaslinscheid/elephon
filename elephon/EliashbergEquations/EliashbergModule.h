/*	This file EliashbergModule.h is part of elephon.
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

#ifndef ELEPHON_ELIASHBERGEQUATIONS_ELIASHBERGMODULE_H_
#define ELEPHON_ELIASHBERGEQUATIONS_ELIASHBERGMODULE_H_

#include "IOMethods/InputOptions.h"
#include <memory>

namespace elephon
{
namespace EliashbergEquations
{

class EliashbergModuleConstants;
class EliashbergFindTc;

/**
 * This is the main driver of the Isotropic Eliashberg equations.
 */
class EliashbergModule
{
public:

	typedef float EliashbergDataType;

	/**
	 * Set up the Eliashberg module dependent on input parameters.
	 *
	 * @param[in] opts
	 */
	EliashbergModule(IOMethods::InputOptions const & opts);

	/**
	 * Dependent on input parameters perform the scheduled duties.
	 */
	void do_work();

	/**
	 * Get the data which remains constant for the entire Eliashberg run schedule
	 * @return a shared pointer to the constant data class.
	 */
	std::shared_ptr<const EliashbergModuleConstants> get_module_constant_data();
private:

	/// for a fixed set of temperatures, run exactly those
	std::vector<double> temperatureSchedule_;

	/// Pointer to a class that handles the determination of Tc if requested by the user.
	std::shared_ptr<EliashbergFindTc> tcFinder_ = nullptr;

	/// Pointer to a class to stores data which remains constant for the entire Eliashberg run schedule.
	std::shared_ptr<EliashbergModuleConstants> moduleConstants_  = nullptr;
};

} /* namespace EliashbergEquations */
} /* namespace elephon */

#endif /* ELEPHON_ELIASHBERGEQUATIONS_ELIASHBERGMODULE_H_ */
