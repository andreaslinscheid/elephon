/*	This file EliashbergModule.cpp is part of elephon.
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

#include "EliashbergEquations/EliashbergModule.h"
#include "EliashbergEquations/EliashbergModuleConstants.h"
#include "PhononStructure/AlphaSquaredF.h"
#include "EliashbergEquations/EliashbergSingleRun.h"
#include "EliashbergEquations/EliashbergSingleRunData.h"
#include "EliashbergEquations/EliashbergFindTc.h"

namespace elephon {
namespace EliashbergEquations {

EliashbergModule::EliashbergModule(IOMethods::InputOptions const & opts)
{
	moduleConstants_ = std::make_shared<EliashbergModuleConstants>();

	// load the a2F iostropic coupling at the Fermi level
	PhononStructure::AlphaSquaredF a2F;
	a2F.load_from_file( opts.get_Eli_f_a2F() );
	moduleConstants_->initialize(opts,std::move(a2F));

	if ( opts.get_EliT().compare("findTc") == 0 ){
		tcFinder_ = std::make_shared<EliashbergFindTc>(opts.get_EliTcThr(), moduleConstants_);
	} else {
		std::stringstream ss(opts.get_EliT());
		double temperature;
		while( ss >> temperature ){
			temperatureSchedule_.push_back(temperature);
		}
	}
}

void
EliashbergModule::do_work()
{
	for (auto t : temperatureSchedule_)
	{
		auto runDataSingle_ptr = std::make_shared<EliashbergSingleRunData>(t, moduleConstants_);
		EliashbergSingleRun run(runDataSingle_ptr);
		bool success = run.solve_Eliashberg();
		if (! success ){
			std::cout << "Failed to converge the Eliashberg equations at temperature " << t << "\n";
		}
	}

	if (tcFinder_)
	{
		tcFinder_->find_Tc();
	}
}

} /* namespace EliashbergEquations */
} /* namespace elephon */
