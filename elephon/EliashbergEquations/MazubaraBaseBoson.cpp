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
		int nMatsBoson,
		int nMatsFermi,
		int numDataPerIndex)
{
	assert(numDataPerIndex > 0);
	nMatsBoson_ = nMatsBoson;
	nMatsFermi_ = nMatsFermi;
	assert(nMatsBoson_ == nMatsFermi_*2+1);
	nDataPerMats_ = numDataPerIndex;
	data_.assign( nMatsBoson_*nDataPerMats_ , T(0.0));
}

void
MatsubaraBaseBoson::set_data(
		int matsubaraIndex,
		int bandIndex,
		int bandPrimeIndex,
		double data)
{
	assert( (matsubaraIndex >= 0) and (nMatsBoson_));
	assert( (bandIndex >= 0) and (bandIndex < nDataPerMats_));
	assert( (bandPrimeIndex >= 0) and (bandPrimeIndex < nDataPerMats_));
	assert( (this->matsubara_to_data_index(matsubaraIndex, bandIndex, bandPrimeIndex) >= 0)
			and (this->matsubara_to_data_index(matsubaraIndex, bandIndex, bandPrimeIndex) < data_.size()));
	data_[this->matsubara_to_data_index(matsubaraIndex, bandIndex, bandPrimeIndex)] = data;
}

void
MatsubaraBaseBoson::get_row_Matsubara_fermi(
		int matsubaraFermiIndex,
		int bandIndex,
		int bandPrimeIndex,
		T const * __restrict & begin,
		T const * __restrict & end) const
{
	// matsubaraFermiIndex = 0 corresponds to the matsubara frequency -nMatsFermi_/2
	// we return the range n' = [-nMatsFermi_/2,nMatsFermi_/2[
	// the data layout is (n'-n)
	begin = data_.data() + (nMatsFermi_ - matsubaraFermiIndex)+nMatsBoson_*(bandPrimeIndex + nDataPerMats_*bandIndex);
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
MatsubaraBaseBoson::matsubara_to_data_index(
		int n,
		int bandIndex,
		int bandPrimeIndex) const
{
	int idata = nDataPerMats_*bandIndex + bandPrimeIndex;
	return idata*nMatsBoson_ + n;
}

} /* namespace EliashbergEquations */
} /* namespace elephon */
