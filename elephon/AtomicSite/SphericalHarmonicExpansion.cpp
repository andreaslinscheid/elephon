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
SphericalHarmonicExpansion::transform(
		symmetry::SymmetryOperation const & sop)
{
	rgrid_.transform(sop);
	sop.rotate_radial_scalar_data(lmax_, rMax_, data_);
}

void
SphericalHarmonicExpansion::set_center(std::vector<double> center)
{
	rgrid_.set_center(std::move(center));
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
