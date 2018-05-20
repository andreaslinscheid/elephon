/*	This file EliashbergFindTc.cpp is part of elephon.
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

#include "EliashbergEquations/EliashbergFindTc.h"
#include "EliashbergEquations/EliashbergModuleConstants.h"
#include "EliashbergEquations/EliashbergSingleRunData.h"
#include "EliashbergEquations/EliashbergSingleRun.h"
#include <vector>
#include <cmath>

namespace elephon {
namespace EliashbergEquations {

EliashbergFindTc::EliashbergFindTc(
		double tcThreshold,
		std::shared_ptr<EliashbergModuleConstants> runConstant_ptr)
	: tcThreshold_(tcThreshold), runConstant_(std::move(runConstant_ptr))
{
}

void
EliashbergFindTc::find_Tc()
{
	// solve McMillan to get a first estimate
	double lambda = runConstant_->get_isotropic_elphon_coupling().compute_coupling_lambda();
	double muStar = *std::max_element(runConstant_->get_isotropic_coulomb().begin(), runConstant_->get_isotropic_coulomb().end());
	double omegaLn = runConstant_->get_isotropic_elphon_coupling().compute_omegaLog(lambda);
	double currentTemperature = tcThreshold_;
	if ( lambda > muStar )
	{
		double McMillanTc = omegaLn / 1.2 * std::exp( - 1.04*(1.0+lambda) / (lambda-muStar*(1+0.62*lambda)) )
							/Auxillary::units::BOLTZMANN_CONSTANT_IN_EV_PER_K / Auxillary::units::EV_TO_THZ_CONVERSION_FACTOR ;

		std::cout << "Starting from half the maximal guess for the critical temperature = " << McMillanTc << '\n';
		currentTemperature = McMillanTc/2.0;
	}
	else
	{
		std::cout << "Guess for the critical temperature is zero, start at almost the lowest possible temperature = " << currentTemperature << '\n';
	}

	auto runData_ptr = std::make_shared<EliashbergSingleRunData>(currentTemperature, runConstant_);
	EliashbergSingleRun startRun(runData_ptr);
	while ( ! startRun.solve_Eliashberg() )
	{ // Failure to converge probably means we are very close to Tc
		currentTemperature /= 2.0;
		if ( currentTemperature < tcThreshold_)
		{
			std::cout << "Hit the minimal temperature threshold:\n";
			std::cout << " Tc = 0.0 K" << "\n";
			return;
		}
		runData_ptr->reset_gap();
		runData_ptr->reset_temperature(currentTemperature);
	}

	// the first attempt is chosen conservative. If we still get a zero gap, let's try the minimal temperature
	// even though this calculation is slow.
	if ( (currentTemperature > tcThreshold_)
			&& (startRun.get_maximal_gap_magnitude() < runConstant_->get_considered_zero_threshold()))
	{
		runData_ptr->reset_gap();
		currentTemperature = tcThreshold_;
		runData_ptr->reset_temperature(currentTemperature);
		startRun.solve_Eliashberg();
	}

	// early exit
	if ( startRun.get_maximal_gap_magnitude() < runConstant_->get_considered_zero_threshold())
	{
		std::cout << "Hit the minimal temperature threshold:\n";
		std::cout << " Tc = 0.0 K" << "\n";
		return;
	}

	// at this point we have a finite resulting gap function. Find tc.
	double currentGuess = currentTemperature;
	for (int niter = 0 ; niter < runConstant_->get_max_number_iterations(); ++niter)
	{
		double nextGuessTemp = this->guess_next_temperature(startRun.get_maximal_gap_magnitude(),currentTemperature);
		if ( std::abs(currentGuess-nextGuessTemp) < tcThreshold_){
			std::cout << "Elishberg equations after "<< niter<<" temperature steps determined to have a Tc of:";
			std::cout << " Tc = " << nextGuessTemp << " K\n";
			return;
		}
		currentTemperature = nextGuessTemp + (currentTemperature-nextGuessTemp)*0.5;
		runData_ptr->reset_temperature(currentTemperature);
		if ( ! startRun.solve_Eliashberg() ){
			std::cout << "Elishberg equations failed to converge. Probably too close to Tc. The current temperature is\n";
			std::cout << " Tc = " << currentTemperature << " K\n";
			return;
		}
		currentGuess = nextGuessTemp;
	}
}

double
EliashbergFindTc::guess_next_temperature(
		T maxGap,
		double temperature)
{
	//add the new data to the set
	temperatureFirstVsMaxDeltaSecond_.push_back(std::make_pair(temperature,maxGap));

	//get the up to two highest pairs below Tc and check if there is a upper bound for Tc
	std::vector< std::pair<double,T> > upToTwoPairsBelowTc, lowestPairAboveTc;
	this->find_up_to_two_largest_T_pairs_smaller_Tc(upToTwoPairsBelowTc);
	this->find_lowest_temperature_pair_above_tc_if_present(lowestPairAboveTc);

	//upToTwoPairsBelowTc empty? then the given data was above Tc and serves as an upper bound
	if ( upToTwoPairsBelowTc.empty() )
	{
		if ( lowestPairAboveTc.empty() )
		{
			throw std::logic_error("Internal error, the pair added did go nowhere.");
		}
		else
		{
			//return max(current temperature by two, minimal allowed temperature)
			double temp = temperature/2;
			if ( temp < tcThreshold_)
				return tcThreshold_;
			return temp;
		}
	}

	//upToTwoPairsBelowTc has only one entry: attempt to choose a temperature
	//		so that we get the second one and are being able to do heuristics
	if ( upToTwoPairsBelowTc.size() == 1 )
	{
		//case that there is no upper bound: just add 1K
		if ( lowestPairAboveTc.empty() )
		{
			return upToTwoPairsBelowTc[0].first + 1.0;
		}
		else
		{
			//in case there is an upper bound, half the interval
			return (upToTwoPairsBelowTc[0].first + lowestPairAboveTc[0].first)/2.0;
		}
	}

	//with two pairs below Tc we perform heuristics
	double newTemperature = this->temperature_heuristic_square_root(upToTwoPairsBelowTc[1],upToTwoPairsBelowTc[0]);

	//if there is no upper bound return
	if ( lowestPairAboveTc.empty() ){
		return newTemperature;
	}

	//if it turns out the temperature guess is above a known non-superconducting temperature half the interval
	if ( newTemperature >  lowestPairAboveTc[0].first - 1e-8 ) {
		return (upToTwoPairsBelowTc[0].first + lowestPairAboveTc[0].first)/2.0;
	}
	return newTemperature;
}

void
EliashbergFindTc::find_up_to_two_largest_T_pairs_smaller_Tc (
		std::vector< std::pair<double,T> > & orderedPairOfPairsWithLargestFirst) const
{
	//find up to two pairs below tc
	std::pair<double,T> highestPairBelowTc = std::make_pair(0.0,T(0.0));
	std::pair<double,T> secondHighestPairBelowTc = std::make_pair(0.0,T(-2.0));
	for ( auto & p : temperatureFirstVsMaxDeltaSecond_ )
	{
		//if the max gap is significantly not zero, consider this pair
		if ( p.second > runConstant_->get_considered_zero_threshold() )
		{
			if ( (highestPairBelowTc.first <= 0) || (p.first  > highestPairBelowTc.first ))
			{
				secondHighestPairBelowTc = highestPairBelowTc;
				highestPairBelowTc = p;
			}
			else if (p.first  > secondHighestPairBelowTc.first ) {
				secondHighestPairBelowTc = p;
			}
		}
	}
	if (highestPairBelowTc.second > 0)
		orderedPairOfPairsWithLargestFirst.push_back(std::move(highestPairBelowTc));
	if (secondHighestPairBelowTc.second > 0)
		orderedPairOfPairsWithLargestFirst.push_back(std::move(secondHighestPairBelowTc));
}

void
EliashbergFindTc::find_lowest_temperature_pair_above_tc_if_present(
		std::vector< std::pair<double,T> > & pair) const
{
	std::pair<double,T> lowestPairAboveTc = std::make_pair(-1.0,0.0);
	for ( auto const & p : temperatureFirstVsMaxDeltaSecond_)
		if ( p.second <= runConstant_->get_considered_zero_threshold() )
			if ( (lowestPairAboveTc.first <= 0) || (p.first  < lowestPairAboveTc.first ))
				lowestPairAboveTc = p;
	if (lowestPairAboveTc.second > 0)
		pair.push_back(lowestPairAboveTc);
}

double
EliashbergFindTc::temperature_heuristic_square_root(
		std::pair<double,T> const & lower, std::pair<double,T> const & higher) const
{
	//the ratio of two evaluated points \Delta(t1) and \Delta(2) of the assumed function \Delta(t)
	auto tcGuess_f = [] (double tc, double R, double t1, double t2, double alpha){
		double x1_2 = std::pow(t1/tc,2);
		double x2_2 = std::pow(t2/tc,2);
		return R - ((1 - x1_2) + alpha*x1_2*std::sqrt(1 - t1/tc))
				/((1 - x2_2) + alpha*x2_2*std::sqrt(1 - t2/tc));
	};

	//the derivative of the ratio of two evaluated points \Delta(t1) and \Delta(2) of the assumed function \Delta(t)
	auto tcGuess_fprime = [] (double tc, double t1, double t2, double alpha){
		double sqrtOf1MinusT1ByTc = std::sqrt(1 - t1/tc);
		double sqrtOf1MinusT2ByTc = std::sqrt(1 - t2/tc);
		double t1_2 = std::pow(t1,2);
		double t2_2 = std::pow(t2,2);
		double tc_2 = std::pow(tc,2);
		return (t1_2*(-4*tc*sqrtOf1MinusT1ByTc + 4*tc*alpha - 5*t1*alpha))
				/(2*tc_2*sqrtOf1MinusT1ByTc*(tc_2 + t2_2*(-1 + sqrtOf1MinusT2ByTc*alpha)))
				+ ( t2_2*(4*tc*(sqrtOf1MinusT2ByTc - alpha) + 5*t2*alpha)*
						(tc_2 + t1_2*(-1 + sqrtOf1MinusT1ByTc*alpha)))
				 /(2*sqrtOf1MinusT2ByTc*std::pow(std::pow(tc,3) + tc*t2_2*(-1 + sqrtOf1MinusT2ByTc*alpha),2));
	};

	double tcGuess=0;
	double R=higher.second/lower.second;
	double tcNew=higher.first+1e-8;
	while(std::fabs(tcGuess -tcNew) > 1e-8) {
		tcGuess = tcNew;
		tcNew -= tcGuess_f(tcGuess,R,higher.first,lower.first,1)/tcGuess_fprime(tcGuess,higher.first,lower.first,1);
	}
	return tcGuess;
}

EliashbergFindTc::T
EliashbergFindTc::predic_gap_magn_T(std::pair<double,T> const & higher, double T, double Tc) const
{
	double const alpha = 1.0;
	double x1_2 = std::pow(higher.first/Tc,2);
	double x2_2 = std::pow(T/Tc,2);
	EliashbergFindTc::T Delta0 = higher.second / ((1.0 - x1_2) + alpha*x1_2*std::sqrt(1 - higher.first/Tc));
	return Delta0*((1 - x2_2) + alpha*x2_2*std::sqrt(1 - T/Tc));
}

} /* namespace EliashbergEquations */
} /* namespace elephon */
