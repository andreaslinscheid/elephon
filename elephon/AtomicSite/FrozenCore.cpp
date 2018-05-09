/*	This file FrozenCore.cpp is part of elephon.
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
 *  Created on: Apr 27, 2018
 *      Author: A. Linscheid
 */

#include "AtomicSite/FrozenCore.h"

namespace elephon
{
namespace AtomicSite
{

void
FrozenCore::initialize(double corePointCharge,
		std::vector<double> electronicFrozenCoreCharge,
		RadialGrid radialGrid)
{
	corePointCharge_ = corePointCharge;
	electronicFrozenCoreCharge_= std::move(electronicFrozenCoreCharge);
	rgrid_ = std::move(radialGrid);
	coreHartreeZDisplacementDataBuffer_.clear();
}

void
FrozenCore::set_hartree_displ(int ir, double radialIntegralCharge) const
{
	// Please note that electronicFrozenCoreCharge_ includes a factor r^2
	coreHartreeZDisplacementDataBuffer_[ir]
			= - Auxillary::units::HARTREE_TO_EV*corePointCharge_
					  /rgrid_.get_radius(ir)/rgrid_.get_radius(ir) * std::sqrt(4.0*M_PI/3.0) // core point charge
			  + Auxillary::units::HARTREE_TO_EV*std::pow(4*M_PI,1.5)/std::sqrt(3.0)*radialIntegralCharge
			  	  	  /rgrid_.get_radius(ir)/rgrid_.get_radius(ir); // frozen electronic core charge
};

void
FrozenCore::add_contribution_prev_two_intervals (int ir, double & integral, std::vector<double> const & y) const
{
	assert(ir >= 2);
	double r = rgrid_.get_radius(ir);

	const double h0 = rgrid_.get_radius(ir-1)-rgrid_.get_radius(ir-2);
	const double h1 = r-rgrid_.get_radius(ir-1);
	const double hsum = h0 + h1;
	const double hprod = h0 * h1;
	const double h0divh1 = h0 / h1;
	integral += hsum/6.0*(y[ir-2]*(2-1.0/h0divh1) + y[ir-1]*hsum*hsum/hprod + y[ir  ]*(2-h0divh1));
};

void
FrozenCore::compute_core_hartree_analytic_z_displacement() const
{
	const int nR = rgrid_.get_num_R();
	coreHartreeZDisplacementDataBuffer_.resize(nR);
	assert(nR > 2);

	/* Note: We have to compute the integral effectively 2 times if we want
	 * an integration rule that spans two intervals, such as a Simpson rule.
	 * Since we need the partial integral from 0 to r, this means we have
	 * even and odd intervals separately, odd starting at 1 where the integral
	 * is simply a trapezoidal rule
	 */
	double radialIntegralChargeEven = rgrid_.get_radius(0)*electronicFrozenCoreCharge_[0];
	this->set_hartree_displ(0, radialIntegralChargeEven);

	double radialIntegralChargeOdd = radialIntegralChargeEven +
			0.5*(rgrid_.get_radius(1)-rgrid_.get_radius(0))*(electronicFrozenCoreCharge_[1]+electronicFrozenCoreCharge_[0]);
	this->set_hartree_displ(1, radialIntegralChargeOdd);

	// set the even and odd intervals
	for (int ir = 2; ir < nR; ++ir)
	{
		double & radInte = ir % 2 == 0 ?  radialIntegralChargeEven : radialIntegralChargeOdd;
		this->add_contribution_prev_two_intervals(ir, radInte, electronicFrozenCoreCharge_);
		this->set_hartree_displ(ir, radInte);
	}
}

double
FrozenCore::total_electronic_charge() const
{
	assert(electronicFrozenCoreCharge_.size()>=2);
	const int nR = rgrid_.get_num_R();
	assert(electronicFrozenCoreCharge_.size() == nR);
	double charge = rgrid_.get_radius(0)*electronicFrozenCoreCharge_[0];
	for (int ir = 2; ir < nR; ir+=2)
		this->add_contribution_prev_two_intervals(ir,charge, electronicFrozenCoreCharge_);

	if ( nR % 2 == 0 )
		charge += 0.5*(rgrid_.get_radius(nR-1)-rgrid_.get_radius(nR-2))*(electronicFrozenCoreCharge_[nR-1]+electronicFrozenCoreCharge_[nR-2]);
	return charge*4*M_PI;
}

RadialGrid const &
FrozenCore::get_radial_grid() const
{
	return rgrid_;
}

} /* namespace AtomicSite */
} /* namespace elephon */
