/*	This file IsotropicElectronPhononCoupling.cpp is part of elephon.
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

#include "EliashbergEquations/IsotropicElectronPhononCoupling.h"
#include <assert.h>

namespace elephon
{
namespace EliashbergEquations
{

void
IsotropicElectronPhononCoupling::initialize(
		std::string const & filename,
		double temperature,
		double energyCutoff)
{
//	std::ifstream
}

void
IsotropicElectronPhononCoupling::initialize(
		std::vector<double> const & frequencies,
		std::vector<double> const & a2F,
		double temperature,
		double energyCutoff)
{
	assert( a2F.size() == frequencies.size());
	assert( a2F.size() >= 2);
	MatsubaraBaseBoson::initialize(temperature, energyCutoff, 1);
	auto m = MatsubaraBaseBoson::MatsubaraFreq(temperature);

	double step = m(1) - m(0);
	for ( int nv = this->min_mats_freq() ; nv <= this->max_mats_freq(); ++nv)
	{
		double mnnp = step*nv*step*nv;
		if ( mnnp == 0 )
			mnnp = 1e-8;

		double integral = 0;
		int nF = static_cast<int>(frequencies.size());
		for ( int o = 1; o < nF-1 ; o++ )
		{
			integral += 2*a2F[o]*frequencies[o]/( frequencies[o]*frequencies[o] + mnnp)
					*(frequencies[o+1]-frequencies[o-1])/2.0;
		}
		integral += 2*a2F[0]*frequencies[0]/( frequencies[0]*frequencies[0] + mnnp)
				*(frequencies[1]-frequencies[0])/2.0;
		integral += 2*a2F[nF-1]*frequencies[nF-1]/( frequencies[nF-1]*frequencies[nF-1] + mnnp)
				*(frequencies[nF-1]-frequencies[nF-2])/2.0;
		this->set_data(nv, 0, integral);
	}
}

} /* namespace EliashbergEquations */
} /* namespace elephon */
