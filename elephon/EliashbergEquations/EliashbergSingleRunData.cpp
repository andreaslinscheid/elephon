/*	This file EliashbergSingleRunData.cpp is part of elephon.
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

#include "EliashbergEquations/EliashbergSingleRunData.h"
#include "PhononStructure/AlphaSquaredF.h"
#include "EliashbergEquations/EliashbergModuleConstants.h"
#include "EliashbergEquations/Eliashberg_helperfunction.hpp"

namespace elephon {
namespace EliashbergEquations {

EliashbergSingleRunData::EliashbergSingleRunData(
		double temperature,
		std::shared_ptr<const EliashbergModuleConstants> runConstants)
{
	const int nMats = fermi_number_matsubara_frequencies(temperature, runConstants->get_Matsubara_cutoff());
	const int nBands = runConstants->get_isotropic_elphon_coupling().get_num_bands();
	EliashbergGapFunction gapDelta;
	gapDelta.initialize(nMats, nBands, 0.01);
	EliashbergFrequencyRenormalization freqRenormZ;
	freqRenormZ.initialize(nMats, nBands, 1+runConstants->get_isotropic_elphon_coupling().compute_coupling_lambda());
	this->initialize(temperature, std::move(gapDelta), std::move(freqRenormZ), runConstants);
}

void EliashbergSingleRunData::initialize(
		double temperature,
		EliashbergGapFunction gapDelta,
		EliashbergFrequencyRenormalization freqRenormZ,
		std::shared_ptr<const EliashbergModuleConstants> runConstants)
{
	temperature_ = temperature;
	runConstants_ = runConstants;
	auto const & a2F = runConstants_->get_isotropic_elphon_coupling();
	effectiveCoupling_.initialize(a2F, temperature, runConstants_->get_Matsubara_cutoff());
	inverseTemperature_ = Algorithms::helperfunctions::inverse_temperature_eV(temperature);
	gapDelta_ = std::move(gapDelta);
	gapDeltaPrev_ = gapDelta_;
	freqRenormZ_ = std::move(freqRenormZ);
	freqRenormZPrev_ = freqRenormZ_;
}

EliashbergGapFunction const &
EliashbergSingleRunData::get_gap_previous_iteration() const
{
	return gapDeltaPrev_;
}

EliashbergGapFunction const &
EliashbergSingleRunData::get_gap() const
{
	return gapDelta_;
}

EliashbergFrequencyRenormalization const &
EliashbergSingleRunData::get_frequencyRenorm_previous_iteration() const
{
	return freqRenormZPrev_;
}

EliashbergFrequencyRenormalization const &
EliashbergSingleRunData::get_frequencyRenorm() const
{
	return freqRenormZ_;
}

EliashbergGapFunction &
EliashbergSingleRunData::access_gap()
{
	return gapDelta_;
}

EliashbergFrequencyRenormalization &
EliashbergSingleRunData::access_frequencyRenorm()
{
	return freqRenormZ_;
}

std::shared_ptr<const EliashbergModuleConstants>
EliashbergSingleRunData::get_run_constants() const
{
	return runConstants_;
}

IsotropicElectronPhononCoupling const &
EliashbergSingleRunData::get_effective_coupling() const
{
	return effectiveCoupling_;
}

double
EliashbergSingleRunData::get_inverse_temp_beta() const
{
	return inverseTemperature_;
}

std::pair<int,int>
EliashbergSingleRunData::get_coulomb_matsubara_cutoff_index() const
{
	double prefactor = M_PI / inverseTemperature_;
	int lastIndexIn = std::floor((runConstants_->get_isotropic_coulomb_cutoff()-prefactor)/(2.0*prefactor));
	return std::make_pair(-lastIndexIn,lastIndexIn+1);
}


int
EliashbergSingleRunData::get_number_Matsubara_freqs_bose() const
{
	double prefactor = M_PI / inverseTemperature_;
	return std::floor(runConstants_->get_Matsubara_cutoff()/(2.0*prefactor));
}

double
EliashbergSingleRunData::get_mixing_parameter() const
{
	return runConstants_->get_mixing_parameter();
}

void
EliashbergSingleRunData::reset_gap()
{
	std::fill_n(gapDelta_.access_data_ptr(), gapDelta_.get_num_mats()*gapDelta_.get_num_bands(), T(1.0));
	std::fill_n(gapDeltaPrev_.access_data_ptr(), gapDeltaPrev_.get_num_mats()*gapDeltaPrev_.get_num_bands(), T(1.0));
}

void
EliashbergSingleRunData::reset_temperature(double temperature)
{
	auto splineMatrix_ptr = std::make_shared<Auxillary::alignedvector::DV>();
	const int nMatsOld = fermi_number_matsubara_frequencies(temperature_, runConstants_->get_Matsubara_cutoff());
	std::vector<double> oldMatsubaraFrequencies(nMatsOld);
	for (int n = 0 ; n < nMatsOld; ++n )
		oldMatsubaraFrequencies[n] = fermi_matsubara_frequency_of_index(n,nMatsOld,inverseTemperature_);

	const int nMatsNew = fermi_number_matsubara_frequencies(temperature, runConstants_->get_Matsubara_cutoff());
	std::vector<double> newMatsubaraFrequencies(nMatsNew);
	for (int n = 0 ; n < nMatsNew; ++n )
		newMatsubaraFrequencies[n] = fermi_matsubara_frequency_of_index(n,nMatsOld,Algorithms::helperfunctions::inverse_temperature_eV(temperature));

	gapDelta_.change_temperature(oldMatsubaraFrequencies, newMatsubaraFrequencies, splineMatrix_ptr);
	freqRenormZ_.change_temperature(oldMatsubaraFrequencies, newMatsubaraFrequencies, splineMatrix_ptr);
	this->initialize(temperature, std::move(gapDelta_), std::move(freqRenormZ_), runConstants_);
}

void
EliashbergSingleRunData::set_next_iteration()
{
	std::swap(gapDelta_, gapDeltaPrev_);
	std::swap(freqRenormZ_, freqRenormZPrev_);
}

double
EliashbergSingleRunData::get_temperature() const
{
	return temperature_;
}


} /* namespace EliashbergEquations */
} /* namespace elephon */
