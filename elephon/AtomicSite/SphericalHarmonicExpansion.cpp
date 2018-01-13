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

namespace elephon
{
namespace AtomicSite
{
void
SphericalHarmonicExpansion::initialize(
		int lmax,
		int rMax,
		Auxillary::alignedvector::ZV data,
		RadialGrid rgrid)
{
	assert(data.size() == rMax*((lmax+1)*(lmax+1)));
	assert(rMax == rgrid.get_num_R());
	lmax_ = lmax;
	rMax_ = rMax;
	data_ = std::move(data);
	rgrid_ = std::move(rgrid);
}

int
SphericalHarmonicExpansion::angular_momentum_layout(
		int l,
		int m) const
{
	assert((m>=-l)&&(m<=l));
	// sum l'=0,l (2l'+1) = (l-1)l + l = l*l
	return l+m + (l*l);
}

std::complex<double> &
SphericalHarmonicExpansion::operator() (int r, int m, int l)
{
	assert((data_.size() > r+rMax_*this->angular_momentum_layout(l,m))
			&& (r+rMax_*this->angular_momentum_layout(l,m) >= 0));
	return data_[r+rMax_*this->angular_momentum_layout(l,m)];
}

std::complex<double>
SphericalHarmonicExpansion::operator() (int r, int m, int l) const
{
	assert((data_.size() > r+rMax_*this->angular_momentum_layout(l,m))
			&& (r+rMax_*this->angular_momentum_layout(l,m) >= 0));
	return data_[r+rMax_*this->angular_momentum_layout(l,m)];
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
	{
		double x = coordinates[ip*3+0];
		double y = coordinates[ip*3+1];
		double z = coordinates[ip*3+2];
		double r = std::sqrt(std::pow(x,2) + std::pow(y,2) + std::pow(z,2));
		spherical_r[ip] = r;
		if (r < 1e-10) // by convention, set theta and phi to zero for zero radius
		{
			spherical_phi[ip] = 0.0;
			spherical_theta[ip] = 0.0;
			continue;
		}
		double theta = std::acos(z/r);
		double phi = std::atan2(y,x);
		spherical_theta[ip] = theta;
		spherical_phi[ip] = phi;
	}

	// for each l and m interpolate the radial data to the r-values and multiply the spherical harmonics
	interpolated_data.assign(numPts, std::complex<double>(0.0));
	auto radial_interpol = interpolated_data;
	for (int l = 0 ; l <= lmax_; ++l)
		for (int m = -l ; m <= l; ++m)
		{
			auto itBegin = data_.begin() + this->angular_momentum_layout(l,m)*rMax_;
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
SphericalHarmonicExpansion::apply_wigner_D_rotation(
		std::vector<WignerDMatrix> const & wignerD)
{
	Algorithms::LinearAlgebraInterface linalg;
	assert(wignerD.size() > lmax_);

	elephon::Auxillary::alignedvector::ZV buffer(rMax_*(2*lmax_+1));
	for (int l = 0 ; l <= lmax_ ; ++l)
	{
		int nl = 2*l+1;
		assert(wignerD[l].view_as_matrix().size() == nl*nl);
		auto wignerD_ptr = wignerD[l].view_as_matrix().data();
		auto expansionData_ptr = data_.data() + this->angular_momentum_layout(l,-l)*rMax_;
		linalg.call_gemm('n', 'n',
				nl, rMax_, nl, // W is a (2*l+1) x (2*l+1) matrix, the expansion data we interpret as a (2*l+1) x rMax_ matrix.
				decltype(*wignerD_ptr)(1.0), wignerD_ptr, nl,
				expansionData_ptr, rMax_,
				decltype(*wignerD_ptr)(0.0), buffer.data(), rMax_ );

		std::copy(buffer.begin(), buffer.begin()+rMax_*nl, expansionData_ptr);
	}
}

} /* namespace AtomicSite */
} /* namespace elephon */
