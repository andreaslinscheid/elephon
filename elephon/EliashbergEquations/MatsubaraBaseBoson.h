/*	This file MatsubaraBaseBoson.h is part of elephon.
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

#ifndef ELEPHON_ELIASHBERGEQUATIONS_MATSUBARABASEBOSON_H_
#define ELEPHON_ELIASHBERGEQUATIONS_MATSUBARABASEBOSON_H_

#include "EliashbergEquations/EliashbergModule.h"
#include "Auxillary/AlignedVector.h"
#include <vector>

namespace elephon
{
namespace EliashbergEquations
{

class MatsubaraBaseBoson
{
	typedef EliashbergModule::EliashbergDataType T;
public:

	void initialize(
			int nMatsBoson,
			int nMatsFermi,
			int numDataPerIndex = 1);

	void set_data(
			int matsubaraIndex,
			int bandIndex,
			int bandPrimeIndex,
			double data);

	void get_row_Matsubara_fermi(
			int matsubaraFermiIndex,
			int bandIndex,
			int bandPrimeIndex,
			T const * __restrict & begin,
			T const * __restrict & end) const;

	int min_mats_freq() const;

	int max_mats_freq() const;

private:

	int nMatsBoson_ = 0;

	int nMatsFermi_ = 0;

	int nDataPerMats_ = 1;

	Auxillary::alignedvector::aligned_vector<T> data_;

	int matsubara_to_data_index(
			int n,
			int bandIndex,
			int bandPrimeIndex) const;
};

namespace detail
{
class MatsubaraFreq
{
public:


};
}; /* namespace detail */

} /* namespace EliashbergEquations */
} /* namespace elephon */

#endif /* ELEPHON_ELIASHBERGEQUATIONS_MATSUBARABASEBOSON_H_ */
