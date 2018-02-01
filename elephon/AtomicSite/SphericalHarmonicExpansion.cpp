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
#include "Algorithms/helperfunctions.hpp"
#include "symmetry/SymmetryOperation.h"

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

std::complex<double> &
SphericalHarmonicExpansion::operator() (int r, int m, int l)
{
	using Auxillary::memlayout::angular_momentum_layout;
	assert((data_.size() > r+rMax_*angular_momentum_layout(l,m))
			&& (r+rMax_*angular_momentum_layout(l,m) >= 0));
	return data_[r+rMax_*angular_momentum_layout(l,m)];
}

std::complex<double>
SphericalHarmonicExpansion::operator() (int r, int m, int l) const
{
	using Auxillary::memlayout::angular_momentum_layout;
	assert((data_.size() > r+rMax_*angular_momentum_layout(l,m))
			&& (r+rMax_*angular_momentum_layout(l,m) >= 0));
	return data_[r+rMax_*angular_momentum_layout(l,m)];
}

void
SphericalHarmonicExpansion::interpolate(
		std::vector<double> const & coordinates,
		Auxillary::alignedvector::ZV & interpolated_data) const
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
	interpolated_data.assign(numPts, std::complex<double>(0.0));
	auto radial_interpol = interpolated_data;
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

} /* namespace AtomicSite */
} /* namespace elephon */
