/*	This file FrozenCore.hpp is part of elephon.
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
#include "Auxillary/UnitConversion.h"
#include "Algorithms/rotation_matrix_connecting_unit_vectors.h"
#include "Auxillary/AlignedVector.h"
#include "AtomicSite/EulerAngles.h"

namespace elephon
{
namespace AtomicSite
{
template<class iterator>
void
FrozenCore::add_core_potential(
		iterator sphericalPotentialToBeAddedToBegin,
		iterator sphericalPotentialToBeAddedToEnd) const
{
	const int nR = std::distance(sphericalPotentialToBeAddedToBegin,sphericalPotentialToBeAddedToEnd);
	assert( nR == rgrid_.get_num_R());
	assert(corePointCharge_ > 0);
	assert(electronicFrozenCoreCharge_.size() == rgrid_.get_num_R());
	for (int ir = 0 ; ir < nR; ++ir)
	{
		// core point charge spherical harmonic expansion coefficients:
		*sphericalPotentialToBeAddedToBegin += -Auxillary::units::HARTREE_TO_EV
				*corePointCharge_/rgrid_.get_radius(ir)*std::sqrt(4*M_PI);
		++sphericalPotentialToBeAddedToBegin;
	}
}

template<class iterator>
void
FrozenCore::add_core_hartree_potential(
		iterator sphericalPotentialToBeAddedToBegin,
		iterator sphericalPotentialToBeAddedToEnd) const
{
	const int nR = std::distance(sphericalPotentialToBeAddedToBegin,sphericalPotentialToBeAddedToEnd);
	assert( nR == rgrid_.get_num_R());
	assert(corePointCharge_ > 0);
	assert(electronicFrozenCoreCharge_.size() == rgrid_.get_num_R());

	// for the frozen electronic core charge potential, we need the partial integrals from 0 to r and from r to infinity for each r
	// we pre-calculate it here. The formula is taken from https://github.com/scipy/scipy/blob/v0.14.0/scipy/integrate/quadrature.py#L315
	auto simpsIntegrator = [&](std::vector<double> &partialInt, std::vector<double> const & y){
		partialInt[0] = rgrid_.get_radius(0)*y[0];

		partialInt[1] = partialInt[0] +
				0.5*(rgrid_.get_radius(1)-rgrid_.get_radius(0))*(y[1]+y[0]);

		// set the even/odd intervals
		for (int ir = 2; ir < nR; ++ir)
		{
			partialInt[ir] = partialInt[ir-2];
			this->add_contribution_prev_two_intervals(ir,partialInt[ir],y);
		}
	};

	// first obtain the partial integral with r^2
	// please note that electronicFrozenCoreCharge_ includes a factor r^2
	std::vector<double> partialInteRSquared(nR);
	simpsIntegrator(partialInteRSquared, electronicFrozenCoreCharge_);

	std::vector<double> cTimesR = electronicFrozenCoreCharge_;
	for (int ir = 0 ; ir < nR; ++ir)
		cTimesR[ir] /= rgrid_.get_radius(ir);
	std::vector<double> partialInteR(nR);
	simpsIntegrator(partialInteR, cTimesR);

	auto pit = sphericalPotentialToBeAddedToBegin;
	for (int ir = 0 ; ir < nR; ++ir)
	{
		const double r = rgrid_.get_radius(ir);
		// frozen electronic core charge Hartree potential spherical harmonic expansion coefficients:
		const double lastmod2 = partialInteR[nR-1];
		const double IntegralRToInfinity = lastmod2 - partialInteR[ir];
		const double Integral0ToR = partialInteRSquared[ir];
		*pit += Auxillary::units::HARTREE_TO_EV*std::pow(4*M_PI,1.5)*(IntegralRToInfinity + Integral0ToR / r);
		++pit;
	}
}

template<class iterator>
void
FrozenCore::add_core_hartree_analytic_displacement(
		std::array<double,3> direction,
		iterator potentialToBeAddedToBegin,
		iterator potentialToBeAddedToEnd) const
{
	const int l = 1; // Angular momentum channel we are taling is 1
	const int nR = rgrid_.get_num_R();
	assert(std::distance(potentialToBeAddedToBegin,potentialToBeAddedToEnd) % nR == 0 );
	assert(electronicFrozenCoreCharge_.size() == rgrid_.get_num_R());
	assert(std::abs(std::pow(direction[0],2)+std::pow(direction[1],2)+std::pow(direction[2],2)-1.0)<1e-3);

	// build the data for the z-direction, then rotate the expansion coefficients to the requested direction
	this->compute_core_hartree_analytic_z_displacement();
	assert(coreHartreeZDisplacementDataBuffer_.size()==nR);

	Auxillary::Multi_array<double, 2> rotMat;
	Algorithms::rotation_matrix_connecting_unit_vectors({0.0,0.0,1.0}, std::vector<double>(direction.begin(), direction.end()), rotMat);
	double alpha,beta,gamma;
	AtomicSite::eulerAngles(std::vector<double>(rotMat.data(), rotMat.data()+rotMat.size()), alpha,beta,gamma);
	WignerDMatrix rotMatSH;
	rotMatSH.initialize(l, alpha, beta, gamma, true);

	// In principle this is a matrix multiplication. Note, however, that coreHartreeZDisplacementDataBuffer_ lives entirely in the
	// l=1,m=0 channel and is real, thus we can simplify this a lot. Note also that the wigner matrix is transposed here
	for (int m = -l; m <=l ;++m)
	{
		auto pit = potentialToBeAddedToBegin + nR*Auxillary::memlayout::angular_momentum_layout(l,m);
		for (int ir = 0; ir < nR; ++ir,++pit)
			*pit += rotMatSH(0,m)*coreHartreeZDisplacementDataBuffer_[ir];
	}
}

} /* namespace AtomicSite */
} /* namespace elephon */
