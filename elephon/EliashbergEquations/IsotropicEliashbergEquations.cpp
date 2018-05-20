/*	This file IsotropicEliashbergEquations.cpp is part of elephon.
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

#include "EliashbergEquations/IsotropicEliashbergEquations.h"
#include "EliashbergEquations/EliashbergSingleRunData.h"
#include "EliashbergEquations/EliashbergModuleConstants.h"
#include "EliashbergEquations/Eliashberg_helperfunction.hpp"
#include <algorithm>

namespace elephon {
namespace EliashbergEquations {

IsotropicEliashbergEquations::IsotropicEliashbergEquations(std::shared_ptr<EliashbergSingleRunData> runData)
{
	runData_ = runData;
}

bool
IsotropicEliashbergEquations::converge(int maximalNumberOfIterations)
{
	int iteration = 0;
	for ( ; iteration < maximalNumberOfIterations; ++iteration)
	{
		// solve the frequency renormalization equation
		this->frequRenorm_equation();

		// then solve the gap equation
		this->gap_equation();

		if ( this->check_convergence() )
		{
			std::cout << "Max(Gap)="<<
					(*std::max_element(runData_->get_gap().get_data_ptr(), runData_->get_gap().get_end_data_ptr()))*1000
					<< "meV at T=" << runData_->get_temperature() << "K"
					" converged in " << (iteration+1) << " iterations.\n";
			return true;
		}

		this->gap_analysis_and_convergence_improvement();

		this->mix_iterations();

		runData_->set_next_iteration();
	}
	return false;
}

void
IsotropicEliashbergEquations::gap_analysis_and_convergence_improvement()
{
	typedef EliashbergModule::EliashbergDataType T;
	T const * __restrict gap_ptr = runData_->get_gap().get_data_ptr();
	T const * __restrict gap_end_ptr = runData_->access_gap().get_end_data_ptr();
	T const * __restrict gapPrev_ptr = runData_->get_gap_previous_iteration().get_data_ptr();


	int numElem = 0;
	auto itg = gap_ptr;
	auto itgP = gapPrev_ptr;
	T averageGap = T(0), averageGapOld = T(0);
	for ( ; itg != gap_end_ptr; ++itg, ++itgP)
	{
		averageGap += *itg;
		averageGapOld += *itgP;
		++numElem;
	}

	averageGap /= numElem;
	averageGapOld /= numElem;
//	std::cout << "Gap analysis:\t Average (new, old) " <<  averageGap << ',' << averageGapOld;

	itg = gap_ptr;
	itgP = gapPrev_ptr;
	T shapeConvergence = 0.0;
	for ( ; itg != gap_end_ptr; ++itg, ++itgP)
		shapeConvergence += std::pow(*itg - (*itgP)*averageGap/averageGapOld, 2);
	shapeConvergence = std::sqrt(shapeConvergence)/ averageGap ;
//	std::cout << "; shape convergence: "<< shapeConvergence <<"\n";

	if ( shapeConvergence < runData_->get_run_constants()->get_shapeConvergence_improvement_threshold())
	{
		T scale = averageGap/averageGapOld;
		T C = runData_->get_run_constants()->get_shapeConvergence_improvement_threshold();
		if ( 100*shapeConvergence < C) {
			for (auto it = runData_->access_gap().access_data_ptr(); it != runData_->access_gap().get_end_data_ptr(); ++it)
				*it *= std::pow(scale,100);
		} else if ( 10*shapeConvergence < C) {
			for (auto it = runData_->access_gap().access_data_ptr(); it != runData_->access_gap().get_end_data_ptr(); ++it)
				*it *= std::pow(scale,10);
		} else if ( shapeConvergence < C) {
			for (auto it = runData_->access_gap().access_data_ptr(); it != runData_->access_gap().get_end_data_ptr(); ++it)
				*it *= std::pow(scale,5);
		}
	}
}

void
IsotropicEliashbergEquations::gap_equation()
{
	typedef EliashbergModule::EliashbergDataType T;
	const int numMatsubaraFermi = runData_->get_gap().get_num_mats();
	const int numBands = runData_->get_run_constants()->get_number_bands();
	const double invBeta = 1.0/runData_->get_inverse_temp_beta();
	const std::pair<int,int> nMatsubaraCoulombCutoff = runData_->get_coulomb_matsubara_cutoff_index();
	std::vector<double> muStar = runData_->get_run_constants()->get_isotropic_coulomb();
	if ((muStar.size() == 1) && (numBands>1))
		muStar.resize(numBands*numBands,muStar[0]);

	EliashbergGapFunction & gap = runData_->access_gap();
	std::fill_n(gap.access_data_ptr(), numBands*numMatsubaraFermi, T(0));

	auto const & freqRenormPrev = runData_->get_frequencyRenorm_previous_iteration();
	auto const & gapPrev = runData_->get_gap_previous_iteration();
	Auxillary::alignedvector::aligned_vector<T> buffer(numBands*numMatsubaraFermi);
	for (int iband = 0 ; iband < numBands; ++iband)
		for (int n = 0 ; n < numMatsubaraFermi; ++n)
			buffer[n + numMatsubaraFermi*iband] = gapPrev(n,iband)/std::sqrt(
					std::pow(fermi_matsubara_frequency_of_index(n,numMatsubaraFermi,runData_->get_inverse_temp_beta()),2)+std::pow(gapPrev(n,iband),2));

	T const * __restrict coupling_ptr_begin = nullptr;
	T const * __restrict coupling_ptr_end = nullptr;
	for (int iband = 0 ; iband < numBands; ++iband)
	{
		for (int ibandPrime = 0 ; ibandPrime < numBands; ++ibandPrime)
			for (int n = 0 ; n < numMatsubaraFermi; ++n)
			{
				// add the phonon coupling
				runData_->get_effective_coupling().get_row_Matsubara_fermi(n,iband,ibandPrime,coupling_ptr_begin, coupling_ptr_end);
				assert(std::distance(coupling_ptr_begin, coupling_ptr_end) == numMatsubaraFermi);
				T const * __restrict buffer_ptr = buffer.data() + ibandPrime*numMatsubaraFermi;

				T sum = T(0);
				for ( int i = 0 ; i < numMatsubaraFermi; ++i, ++coupling_ptr_begin, ++buffer_ptr)
					sum += (*coupling_ptr_begin) * (*buffer_ptr);

				// add mu* repulsion
				const T muStarVal = muStar[iband + numBands*ibandPrime];
				if ( muStarVal != T(0.0) )
				{
					buffer_ptr = buffer.data() + ibandPrime*numMatsubaraFermi;
					for ( int i = nMatsubaraCoulombCutoff.first ; i < nMatsubaraCoulombCutoff.second; ++i, ++coupling_ptr_begin, ++buffer_ptr)
						sum -=  muStarVal * (*buffer_ptr);
				}

				gap(n,iband) += M_PI*sum* invBeta / freqRenormPrev(n,iband);
			}
	}
}

void
IsotropicEliashbergEquations::frequRenorm_equation()
{
	typedef EliashbergModule::EliashbergDataType T;
	const int numMatsubaraFermi = runData_->get_gap().get_num_mats();
	const int numBands = runData_->get_run_constants()->get_number_bands();
	const double invBeta = 1.0/runData_->get_inverse_temp_beta();

	auto & frequencyRenormalization = runData_->access_frequencyRenorm();
	std::fill_n(frequencyRenormalization.access_data_ptr(), numBands*numMatsubaraFermi, T(1.0));
	auto const & gapPrev = runData_->get_gap_previous_iteration();

	Auxillary::alignedvector::aligned_vector<T> buffer(numBands*numMatsubaraFermi);
	for (int iband = 0 ; iband < numBands; ++iband)
		for (int n = 0 ; n < numMatsubaraFermi; ++n)
		{
			auto omegan = fermi_matsubara_frequency_of_index(n,numMatsubaraFermi,runData_->get_inverse_temp_beta());
			buffer[n + numMatsubaraFermi*iband] = omegan / std::sqrt(std::pow(omegan,2)+std::pow(gapPrev(n,iband),2));
		}

	T const * __restrict coupling_ptr_begin = nullptr;
	T const * __restrict coupling_ptr_end = nullptr;
	for (int iband = 0 ; iband < numBands; ++iband)
	{
		for (int ibandPrime = 0 ; ibandPrime < numBands; ++ibandPrime)
			for (int n = 0 ; n < numMatsubaraFermi; ++n)
			{
				// add the phonon coupling
				runData_->get_effective_coupling().get_row_Matsubara_fermi(n,iband,ibandPrime,coupling_ptr_begin, coupling_ptr_end);
				assert(std::distance(coupling_ptr_begin, coupling_ptr_end) == numMatsubaraFermi);
				T const * __restrict buffer_ptr = buffer.data() + ibandPrime*numMatsubaraFermi;

				T sum = T(0);
				for ( int i = 0 ; i < numMatsubaraFermi; ++i, ++coupling_ptr_begin, ++buffer_ptr)
					sum += (*coupling_ptr_begin) * (*buffer_ptr);

				auto omegan = fermi_matsubara_frequency_of_index(n,numMatsubaraFermi,runData_->get_inverse_temp_beta());
				frequencyRenormalization(n,iband) += M_PI*sum*invBeta / omegan;
			}
	}
}

bool
IsotropicEliashbergEquations::check_convergence() const
{
	typedef EliashbergModule::EliashbergDataType T;
	double cth = runData_->get_run_constants()->get_convergence_threshold();
	double zeroEquiv = runData_->get_run_constants()->get_considered_zero_threshold();

	auto cmp_values = [&] (T v1, T v2){
		auto diff = std::max(std::abs(v1),std::abs(v2))*cth;
		if (std::abs(v1-v2)<diff)
			return true;
		if (std::abs(v1)<zeroEquiv && (std::abs(v2)<zeroEquiv))
			return true;
		return false;
	};

	T const * __restrict gap_ptr = runData_->get_gap().get_data_ptr();
	T const * __restrict gap_end_ptr = runData_->get_gap().get_end_data_ptr();
	T const * __restrict gapPrev_ptr = runData_->get_gap_previous_iteration().get_data_ptr();
	for ( ; gap_ptr != gap_end_ptr; ++gap_ptr, ++gapPrev_ptr)
		if (! cmp_values(*gap_ptr,*gapPrev_ptr) )
			return false;

	T const * __restrict Z_ptr = runData_->get_frequencyRenorm().get_data_ptr();
	T const * __restrict Z_end_ptr = runData_->access_frequencyRenorm().get_end_data_ptr();
	T const * __restrict ZPrev_ptr = runData_->get_frequencyRenorm_previous_iteration().get_data_ptr();
	for ( ; Z_ptr != Z_end_ptr; ++Z_ptr, ++ZPrev_ptr)
		if (! cmp_values(*Z_ptr,*ZPrev_ptr) )
			return false;

	return true;
}

void
IsotropicEliashbergEquations::mix_iterations()
{
	typedef EliashbergModule::EliashbergDataType T;
	auto a = runData_->get_mixing_parameter();

	T * __restrict gap_ptr = runData_->access_gap().access_data_ptr();
	T const * __restrict gap_end_ptr = runData_->access_gap().get_end_data_ptr();
	T const * __restrict gapPrev_ptr = runData_->get_gap_previous_iteration().get_data_ptr();

	auto itg = gap_ptr;
	auto itgP = gapPrev_ptr;
	for ( ; itg != gap_end_ptr; ++itg, ++itgP)
		*itg += a*(*itgP - *itg);

	T * __restrict Z_ptr = runData_->access_frequencyRenorm().access_data_ptr();
	T const * __restrict Z_end_ptr = runData_->access_frequencyRenorm().get_end_data_ptr();
	T const * __restrict ZPrev_ptr = runData_->get_frequencyRenorm_previous_iteration().get_data_ptr();

	itg = Z_ptr;
	itgP = ZPrev_ptr;
	for ( ; itg != Z_end_ptr; ++itg, ++itgP)
		*itg += a*(*itgP - *itg);
}

} /* namespace EliashbergEquations */
} /* namespace elephon */
