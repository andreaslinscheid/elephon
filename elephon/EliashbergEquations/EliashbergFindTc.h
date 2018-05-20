/*	This file EliashbergFindTc.h is part of elephon.
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
 *  Created on: May 17, 2018
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_ELIASHBERGEQUATIONS_ELIASHBERGFINDTC_H_
#define ELEPHON_ELIASHBERGEQUATIONS_ELIASHBERGFINDTC_H_

#include "EliashbergEquations/EliashbergModule.h"
#include <memory>
#include <vector>

namespace elephon {
namespace EliashbergEquations {

class EliashbergModuleConstants;
class EliashbergSingleRun;

class EliashbergFindTc
{
	typedef EliashbergModule::EliashbergDataType T;
public:

	EliashbergFindTc( double tcThreshold,
			std::shared_ptr<EliashbergModuleConstants> runConstant_ptr);

	void find_Tc();

private:

	double tcThreshold_;

	std::shared_ptr<EliashbergModuleConstants> runConstant_;

	std::vector<std::pair<double,T>> temperatureFirstVsMaxDeltaSecond_;

	double guess_next_temperature(T maxGap, double currentTemperature);

	void find_up_to_two_largest_T_pairs_smaller_Tc(std::vector< std::pair<double,T> > & pairs) const;

	void find_lowest_temperature_pair_above_tc_if_present(std::vector< std::pair<double,T> > & pairs) const;

	/**
	 * 	Guess a new critical temperature.
	 *
	 * 	The guess is based on the formula \f$\Delta(T) = \Delta_0 [1-(T/T_c)^2+a (T/T_c)^2\sqrt{1-T/T_c}]\f$
	 * 	A new guess requires two points that lead to \f$ T_c(\Delta(T_1)/\Delta(T_2),T_1,T_2)\f$
	 * 	\f$ T_c\f$ has to be found numerically
	 *
	 * @param lower
	 * @param higher
	 * @return
	 */
	double temperature_heuristic_square_root(std::pair<double,T> const & lower, std::pair<double,T> const & higher) const;

	/**
	 * 	Guess a new gap at given temperature.
	 *
	 * 	The guess is based on the formula \f$\Delta(T) = \Delta_0 [1-(T/T_c)^2+a (T/T_c)^2\sqrt{1-T/T_c}]\f$
	 *
	 * @param higher
	 * @param T
	 * @param Tc
	 * @return
	 */
	T predic_gap_magn_T(std::pair<double,T> const & higher, double T, double Tc) const;
};

} /* namespace EliashbergEquations */
} /* namespace elephon */

#endif /* ELEPHON_ELIASHBERGEQUATIONS_ELIASHBERGFINDTC_H_ */
