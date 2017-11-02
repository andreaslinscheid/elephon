/*	This file MatsubaraBaseBoson.cpp is part of elephon.
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
 *  Created on: Oct 25, 2017
 *      Author: A. Linscheid
 */

#include "Algorithms/helperfunctions.hpp"
#include "EliashbergEquations/MatsubaraBaseBoson.h"
#include <assert.h>
#include <math.h>

namespace elephon
{
namespace EliashbergEquations
{

void
MatsubaraBaseBoson::initialize(
		double temperature,
		double energyCutoffMats,
		int numDataPerIndex)
{
	assert(numDataPerIndex > 0);
	nMatsBoson_ = this->compute_n_mats_cutoff(temperature, energyCutoffMats);
//	nMatsFermi_ = Algorithms::helperfunctions::num_Mats_Fermi_cutoff(temperature, energyCutoffMats);
	nDataPerMats_ = numDataPerIndex;
	data_.assign( nMatsBoson_*nDataPerMats_ , 0.0);
}

void
MatsubaraBaseBoson::set_data(
		int mazubaraIndex,
		int dataIndex,
		double data)
{
	assert( (mazubaraIndex >= this->min_mats_freq()) and (mazubaraIndex <= this->max_mats_freq()));
	assert( (dataIndex >= 0) and (dataIndex < nDataPerMats_));
	assert( (this->matsubara_to_data_index(mazubaraIndex, dataIndex) >= 0)
			and (this->matsubara_to_data_index(mazubaraIndex, dataIndex) < data_.size()));
	data_[this->matsubara_to_data_index(mazubaraIndex, dataIndex)] = data;
}

void
MatsubaraBaseBoson::get_row_Matsubara_fermi(
		int mazubaraIndex,
		int dataIndex,
		std::vector<double>::const_iterator & begin,
		std::vector<double>::const_iterator & end) const
{
	begin = data_.begin() + this->matsubara_to_data_index(mazubaraIndex, dataIndex);
	end = begin + nMatsFermi_;
}

int
MatsubaraBaseBoson::min_mats_freq() const
{
	return -nMatsBoson_/2;
}

int
MatsubaraBaseBoson::max_mats_freq() const
{
	return nMatsBoson_/2;
}

int
MatsubaraBaseBoson::compute_n_mats_cutoff(
		double temperature,
		double energyCutoffMats ) const
{
//	return Algorithms::helperfunctions::num_Mats_Fermi_cutoff(temperature, energyCutoffMats)*2 + 1;
	return 0;
}

int
MatsubaraBaseBoson::matsubara_to_data_index(int n, int idata) const
{
	return idata*nMatsBoson_ + n + nMatsBoson_/2;
}

namespace detail
{

MatsubaraFreq::MatsubaraFreq(double temperature)
	: inverseTemperature_( 0.0)// Algorithms::helperfunctions::inverse_temperature_eV(temperature) )
{
}

double
MatsubaraFreq::operator() (int index) const
{
	return 2.0 * M_PI / inverseTemperature_ * index;
}

std::vector<double>
MatsubaraFreq::range(int istart, int iend) const
{
	if ( iend <= istart )
		return std::vector<double>();

	std::vector<double> result;
	result.reserve(iend-istart);
	for ( int n = istart; n < iend ; ++n)
		result.push_back( (*this)(n) );
	return result;

}

} /* namespace detail */
} /* namespace EliashbergEquations */
} /* namespace elephon */
