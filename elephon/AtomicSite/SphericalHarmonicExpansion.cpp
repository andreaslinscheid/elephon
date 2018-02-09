/*	This file SphericalHarmonicExpansion.cpp is part of elephon.
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
 *  Created on: Jan 4, 2018
 *      Author: A. Linscheid
 */

#include "AtomicSite/SphericalHarmonicExpansion.h"
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include "Auxillary/memory_layout_functions.hpp"
#include "symmetry/SymmetryOperation.h"
#include "Algorithms/helperfunctions.hpp"

namespace elephon
{
namespace AtomicSite
{
void
SphericalHarmonicExpansion::initialize(
		int lmax,
		Auxillary::alignedvector::ZV data,
		RadialGrid rgrid)
{
	assert(data.size() == rgrid.get_num_R()*((lmax+1)*(lmax+1)));
	rMax_ = rgrid.get_num_R();
	lmax_ = lmax;
	data_ = std::move(data);
	rgrid_ = std::move(rgrid);
}

void
SphericalHarmonicExpansion::interpolate(
		std::vector<double> const & coordinates,
		Auxillary::alignedvector::ZV & interpolated_data) const
{
	assert(coordinates.size()%3 == 0);
	int numPts = coordinates.size()/3;
	interpolated_data.resize(numPts);
	this->interpolate(coordinates, interpolated_data.data());
}

void
SphericalHarmonicExpansion::interpolate(
		std::vector<double> const & coordinates,
		std::complex<double> * interpolated_data) const
{
	assert(coordinates.size()%3 == 0);

	// calculate the spherical coordinates
	int numPts = coordinates.size()/3;
	std::vector<double> spherical_r(numPts);
	std::vector<double> spherical_phi(numPts);
	std::vector<double> spherical_theta(numPts);
	for (int ip = 0 ; ip < numPts; ++ip)
		Algorithms::helperfunctions::compute_spherical_coords(
				coordinates[ip*3+0], coordinates[ip*3+1], coordinates[ip*3+2],
				spherical_r[ip], spherical_theta[ip], spherical_phi[ip]);

	// for each l and m interpolate the radial data to the r-values and multiply the spherical harmonics
	Auxillary::alignedvector::ZV radial_interpol(numPts, std::complex<double>(0.0));
	for (int l = 0 ; l <= lmax_; ++l)
		for (int m = -l ; m <= l; ++m)
		{
			auto itBegin = data_.begin() + Auxillary::memlayout::angular_momentum_layout(l,m)*rMax_;
			auto itEnd = itBegin + rMax_;
			rgrid_.interpolate(
					spherical_r,
					1,
					itBegin,
					itEnd,
					radial_interpol.begin());

			for (int ip = 0 ; ip < numPts; ++ip)
			{
				auto ylm = boost::math::spherical_harmonic(l, m, spherical_theta[ip], spherical_phi[ip] );
				interpolated_data[ip] += ylm*radial_interpol[ip];
			}
		}
}

void
SphericalHarmonicExpansion::transform(
		symmetry::SymmetryOperation const & sop)
{
	rgrid_.transform(sop);
	sop.rotate_radial_data(lmax_, rMax_, data_);
}

void
SphericalHarmonicExpansion::set_ylm_conjg_data(
		int l, int m,
		Auxillary::alignedvector::DV const & thetas,
		Auxillary::alignedvector::DV const & phis,
		Auxillary::alignedvector::ZV & ylmData) const
{
	assert(thetas.size() == phis.size());
	ylmData.resize(thetas.size());
	for (int ip = 0 ; ip < static_cast<int>(thetas.size()); ++ip)
		ylmData[ip] = std::conj(boost::math::spherical_harmonic(l, m, thetas[ip], phis[ip] ));

}

RadialGrid const &
SphericalHarmonicExpansion::get_radial_grid() const
{
	return rgrid_;
}

int
SphericalHarmonicExpansion::get_l_max() const
{
	return lmax_;
}

} /* namespace AtomicSite */
} /* namespace elephon */
