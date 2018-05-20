/*	This file MatsubaraBaseFermion.cpp is part of elephon.
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
 *  Created on: May 16, 2018
 *      Author: A. Linscheid
 */

#include "EliashbergEquations/MatsubaraBaseFermion.h"
#include "Algorithms/CubeSplineInterpolation.h"

namespace elephon {
namespace EliashbergEquations {

MatsubaraBaseFermion::T *
MatsubaraBaseFermion::access_data_ptr()
{
	return data_.data();
}

MatsubaraBaseFermion::T const *
MatsubaraBaseFermion::get_data_ptr() const
{
	return data_.data();
}

void
MatsubaraBaseFermion::initialize(
		int numMatsubara,
		int numBands,
		Auxillary::alignedvector::aligned_vector<T> data)
{
	nMats_ = numMatsubara;
	nDataPerMats_ = numBands;
	data_ = std::move(data);
}

typename MatsubaraBaseFermion::T &
MatsubaraBaseFermion::operator() (int n, int b)
{
	assert((n>=0)&&(n<nMats_));
	assert((b>=0)&&(b<nDataPerMats_));
	return data_[n + nMats_*b];
}

typename MatsubaraBaseFermion::T const &
MatsubaraBaseFermion::operator() (int n, int b) const
{
	assert((n>=0)&&(n<nMats_));
	assert((b>=0)&&(b<nDataPerMats_));
	return data_[n + nMats_*b];
}

typename MatsubaraBaseFermion::T const *
MatsubaraBaseFermion::get_end_data_ptr() const
{
	return data_.data()+data_.size();
}

int
MatsubaraBaseFermion::get_num_bands() const
{
	return nDataPerMats_;
}

int
MatsubaraBaseFermion::get_num_mats() const
{
	return nMats_;
}

void
MatsubaraBaseFermion::change_temperature(
		std::vector<double> const & oldMatsubaraFrequencies,
		std::vector<double> const & newMatsubaraFrequencies,
		std::shared_ptr<Auxillary::alignedvector::aligned_vector<double>> splineMatrixPtr)
{
	Algorithms::CubeSplineInterpolation<double,double> cspline;
	decltype(data_) buffer( this->get_num_bands() * newMatsubaraFrequencies.size() );
	for (int b = 0 ; b < this->get_num_bands(); ++b)
	{
		cspline.initialize(
				data_.begin() + b*nMats_, data_.begin() + (b+1)*nMats_,
				oldMatsubaraFrequencies.begin(), oldMatsubaraFrequencies.end(),
				oldMatsubaraFrequencies.front(), oldMatsubaraFrequencies.back(),
				splineMatrixPtr);
		for (int iMatsFreq = 0 ; iMatsFreq < newMatsubaraFrequencies.size(); ++iMatsFreq )
			buffer[iMatsFreq + newMatsubaraFrequencies.size()*b] = cspline(newMatsubaraFrequencies[iMatsFreq]);
	}
	this->initialize(newMatsubaraFrequencies.size(), nDataPerMats_, std::move(buffer));
}

} /* namespace EliashbergEquations */
} /* namespace elephon */
