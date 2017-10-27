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

#include <vector>

namespace elephon
{
namespace EliashbergEquations
{

//forward declare
namespace detail{ class MatsubaraFreq; };

class MatsubaraBaseBoson
{
public:

	typedef detail::MatsubaraFreq MatsubaraFreq;

	void initialize(
			double temperature,
			double energyCutoffMats,
			int numDataPerIndex = 1);

	void set_data(
			int mazubaraIndex,
			int dataIndex,
			double data);

	void get_row_Matsubara_fermi(
			int mazubaraIndex,
			int dataIndex,
			std::vector<double>::const_iterator & begin,
			std::vector<double>::const_iterator & end) const;

	int min_mats_freq() const;

	int max_mats_freq() const;

private:

	int nMatsBoson_ = 0;

	int nMatsFermi_ = 0;

	int nDataPerMats_ = 1;

	std::vector<double> data_;

	int compute_n_mats_cutoff(
			double temperature,
			double energyCutoffMats ) const;

	int matsubara_to_data_index(int n, int idata) const;
};

namespace detail
{
class MatsubaraFreq
{
public:

	MatsubaraFreq(double temperature);

	double operator() (int index) const;

	std::vector<double> range(int istart, int iend) const;

private:

	double inverseTemperature_;
};
}; /* namespace detail */

} /* namespace EliashbergEquations */
} /* namespace elephon */

#endif /* ELEPHON_ELIASHBERGEQUATIONS_MATSUBARABASEBOSON_H_ */
