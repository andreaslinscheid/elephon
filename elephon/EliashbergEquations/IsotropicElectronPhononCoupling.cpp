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
#include "Algorithms/SimpsonIntegrator.h"
#include "EliashbergEquations/Eliashberg_helperfunction.hpp"
#include <assert.h>

namespace elephon
{
namespace EliashbergEquations
{

void
IsotropicElectronPhononCoupling::initialize(
		PhononStructure::AlphaSquaredF const & a2F,
		double temperature,
		double energyCutoff)
{
	const int nB = a2F.get_num_bands();
	const int nF = a2F.get_num_frequency_samples();
	const int nMats = fermi_number_matsubara_frequencies(temperature, energyCutoff);
	MatsubaraBaseBoson::initialize(nMats*2+1 , nMats, nB);

	auto const & omega = a2F.get_frequencies();
	Algorithms::SimpsonIntegrator<std::remove_reference<decltype(omega)>::type> si;
	Auxillary::alignedvector::DV buffer(nF);

	double step = 2.0 * M_PI / Algorithms::helperfunctions::inverse_temperature_eV(temperature);
	for ( int nv = this->min_mats_freq() ; nv <= this->max_mats_freq(); ++nv)
	{
		const int MatsubaraIndex = nv - this->min_mats_freq();
		double mnnp = step*nv*step*nv;
		if ( mnnp == 0 )
			mnnp = 1e-8;

		for (int iband = 0 ;iband < nB; ++iband )
			for (int ibandPrime = 0 ;ibandPrime < nB; ++ibandPrime )
			{
				for ( int ifreq = 0 ; ifreq < nF; ++ifreq )
				{
					buffer[ifreq] = 2.0 * a2F(ifreq, iband, ibandPrime) * omega[ifreq]
										/ (std::pow(omega[ifreq],2)+mnnp*std::pow(Auxillary::units::EV_TO_THZ_CONVERSION_FACTOR,2));
				}
				this->set_data(MatsubaraIndex, iband, ibandPrime, si.integrate(buffer.begin(), buffer.end(), omega.begin(), omega.end()));
			}
	}
}

} /* namespace EliashbergEquations */
} /* namespace elephon */
