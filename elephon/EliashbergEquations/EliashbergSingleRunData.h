/*	This file EliashbergSingleRunData.h is part of elephon.
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

#ifndef ELEPHON_ELIASHBERGEQUATIONS_ELIASHBERGSINGLERUNDATA_H_
#define ELEPHON_ELIASHBERGEQUATIONS_ELIASHBERGSINGLERUNDATA_H_

#include "EliashbergEquations/EliashbergGapFunction.h"
#include "EliashbergEquations/EliashbergFrequencyRenormalization.h"
#include "EliashbergEquations/IsotropicElectronPhononCoupling.h"

namespace elephon {
namespace EliashbergEquations {

/**
 * A class that stores the data required for a single run of the isotropic Eliashberg equations
 */
class EliashbergSingleRunData
{
	typedef EliashbergModule::EliashbergDataType T;
public:
	EliashbergSingleRunData(
			double temperature,
			std::shared_ptr<const EliashbergModuleConstants> runConstants);

	void initialize(
			double temperature,
			EliashbergGapFunction gapDelta,
			EliashbergFrequencyRenormalization freqRenormZ,
			std::shared_ptr<const EliashbergModuleConstants> runConstants);

	EliashbergGapFunction const & get_gap_previous_iteration() const;

	EliashbergGapFunction const & get_gap() const;

	EliashbergFrequencyRenormalization const & get_frequencyRenorm_previous_iteration() const;

	EliashbergFrequencyRenormalization const & get_frequencyRenorm() const;

	EliashbergGapFunction & access_gap();

	EliashbergFrequencyRenormalization & access_frequencyRenorm();

	std::shared_ptr<const EliashbergModuleConstants> get_run_constants() const;

	IsotropicElectronPhononCoupling const & get_effective_coupling() const;

	double get_inverse_temp_beta() const;

	/**
	 * Get the first index that is inside- and the first index that is outside of the screened Coulomb cutoff.
	 * @return
	 */
	std::pair<int,int> get_coulomb_matsubara_cutoff_index() const;

	int get_number_Matsubara_freqs_bose() const;

	double get_mixing_parameter() const;

	void reset_gap();

	void reset_temperature(double temperature);

	void set_next_iteration();

	double get_temperature() const;
private:
	std::shared_ptr<const EliashbergModuleConstants> runConstants_;

	EliashbergGapFunction gapDelta_, gapDeltaPrev_;

	EliashbergFrequencyRenormalization freqRenormZ_, freqRenormZPrev_;

	IsotropicElectronPhononCoupling effectiveCoupling_;

	double temperature_ = 0;

	double inverseTemperature_ = 0;
};

} /* namespace EliashbergEquations */
} /* namespace elephon */

#endif /* ELEPHON_ELIASHBERGEQUATIONS_ELIASHBERGSINGLERUNDATA_H_ */
