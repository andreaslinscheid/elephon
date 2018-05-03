/*	This file SphericalHarmonicExpansion.hpp is part of elephon.
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
 *  Created on: Feb 1, 2018
 *      Author: A. Linscheid
 */

#include "AtomicSite/SphericalHarmonicExpansion.h"
#include "Algorithms/helperfunctions.hpp"
#include "Algorithms/SphereIntegrator.h"
#include <boost/math/special_functions/spherical_harmonic.hpp>

namespace elephon
{
namespace AtomicSite
{

template<class VT>
void
SphericalHarmonicExpansion::initialize_from_real(
		int lmax,
		VT const & dataExpansionRealSH,
		RadialGrid rgrid)
{
	using Auxillary::memlayout::angular_momentum_layout;
	typedef Auxillary::alignedvector::ZV::value_type CT;
	const int nR = rgrid.get_num_R();
	assert(std::pow(lmax+1,2)*nR == dataExpansionRealSH.size());

	Auxillary::alignedvector::ZV dataComplex(dataExpansionRealSH.size(), CT(0.0));
	for (int l = 0 ; l < lmax; ++l)
	{
		// m = 0
		for (int ir = 0 ; ir < nR ; ++ir)
			dataComplex[ir + nR*angular_momentum_layout(l,/* m = */0)] =
					dataExpansionRealSH[ir + nR*angular_momentum_layout(l,/* m = */0)];
		// m < 0
		for (int m = -l; m < 0; ++m)
			for (int ir = 0 ; ir < nR ; ++ir)
				dataComplex[ir + nR*angular_momentum_layout(l,m)] = 1.0/std::sqrt(2.0)
						*(CT(dataExpansionRealSH[ir + nR*angular_momentum_layout(l, - m)])+
								CT(0,1.0)*CT(dataExpansionRealSH[ir + nR*angular_momentum_layout(l,m)]));
		// m > 0
		for (int m = 1; m <= l; ++m)
			for (int ir = 0 ; ir < nR ; ++ir)
				dataComplex[ir + nR*angular_momentum_layout(l,m)] =	(m%2==0?1.0:-1.0)/std::sqrt(2.0)
						*(CT(dataExpansionRealSH[ir + nR*angular_momentum_layout(l,m)])
								-CT(0,1.0)*CT(dataExpansionRealSH[ir + nR*angular_momentum_layout(l,-m)]));
	}

	this->initialize(lmax,std::move(dataComplex),std::move(rgrid));
}

inline std::complex<double> &
SphericalHarmonicExpansion::operator() (int r, int m, int l)
{
	using Auxillary::memlayout::angular_momentum_layout;
	assert((data_.size() > r+rMax_*angular_momentum_layout(l,m))
			&& (r+rMax_*angular_momentum_layout(l,m) >= 0));
	return data_[r+rMax_*angular_momentum_layout(l,m)];
}

inline std::complex<double>
SphericalHarmonicExpansion::operator() (int r, int m, int l) const
{
	using Auxillary::memlayout::angular_momentum_layout;
	assert((data_.size() > r+rMax_*angular_momentum_layout(l,m))
			&& (r+rMax_*angular_momentum_layout(l,m) >= 0));
	return data_[r+rMax_*angular_momentum_layout(l,m)];
}

template<class VT>
void
SphericalHarmonicExpansion::interpolate(
		VT const & coordinates,
		Auxillary::alignedvector::ZV & interpolated_data) const
{
	assert(coordinates.size()%3 == 0);
	int numPts = coordinates.size()/3;
	interpolated_data.resize(numPts);
	this->interpolate(coordinates, interpolated_data.data());
}

template<class VT>
void
SphericalHarmonicExpansion::interpolate(
		VT const & coordinates,
		std::complex<double> * interpolated_data) const
{
	assert(coordinates.size()%3 == 0);

	// calculate the spherical coordinates
	int numPts = coordinates.size()/3;
	std::vector<double> spherical_r(numPts);
	std::vector<double> spherical_phi(numPts);
	std::vector<double> spherical_theta(numPts);
	std::vector<double> const & c = rgrid_.get_center();
	for (int ip = 0 ; ip < numPts; ++ip)
		Algorithms::helperfunctions::compute_spherical_coords(
				coordinates[ip*3+0] - c[0], coordinates[ip*3+1]- c[1], coordinates[ip*3+2] - c[2],
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

template<class F>
void
SphericalHarmonicExpansion::fit_to_data(
		F const & f,
		int lmax,
		RadialGrid rgrid )
{
	// this is the functor that is integrated on a grid on the unit sphere.
	struct SphericalPart
	{
		SphericalPart(
				F const & f,
				double r0x, double r0y, double r0z,
				Auxillary::alignedvector::ZV const & ylmData)
				: f_(f), r0x_(r0x),r0y_(r0y),r0z_(r0z),ylmData_(ylmData){};

		void set_radius(double radius) {r_ = radius;};

		 // only used to match the signature
		std::complex<double> operator()(double theta, double phi) const;

		// this is the real function used to produce the data on the sphere for given radius.
		void evaluate_many(
				double const * thetas,
				double const * phis,
				int numEvals,
				std::complex<double> * fvalues) const
		{
			assert(numEvals == ylmData_.size());
			coordsBuffer_.resize(numEvals*3);
			for (int iPt = 0; iPt < numEvals; ++iPt)
			{
				coordsBuffer_[iPt*3+0] = r_* std::sin(thetas[iPt])*std::cos(phis[iPt]) + r0x_;
				coordsBuffer_[iPt*3+1] = r_* std::sin(thetas[iPt])*std::sin(phis[iPt]) + r0y_;
				coordsBuffer_[iPt*3+2] = r_* std::cos(thetas[iPt])                     + r0z_;
			}
			f_.interpolate(coordsBuffer_, fvalues);

			// multiply by the ylm data
			for (int iPt = 0; iPt < numEvals; ++iPt)
				fvalues[iPt] *= ylmData_[iPt];
		}

		F const & f_;

		double r_ = 0.0;

		const double r0x_, r0y_, r0z_;

		Auxillary::alignedvector::ZV const & ylmData_;

		mutable Auxillary::alignedvector::DV coordsBuffer_;
	}; // end struct SphericalPart

	Algorithms::SphereIntegrator<SphericalPart> integrator;
	auto lebedev_rule = integrator.pick_rule_spherical_harmonic(2*lmax);
	integrator.initialize(lebedev_rule);
	const int nR = rgrid.get_num_R();
	const int nSP = integrator.get_num_pts_surface();

	Auxillary::alignedvector::DV phis, thetas;
	integrator.get_surface_pts(thetas, phis);

	Auxillary::alignedvector::ZV coefficientData((lmax+1)*(lmax+1)*nR);
	Auxillary::alignedvector::ZV ylmData(nSP);

	auto c = rgrid.get_center();
	for (int iL = 0 ; iL < lmax; ++iL )
	{
		for (int iM = -iL ; iM <= iL; ++iM )
		{
			// we pre-compute the YLM data since it is the same for all
			// radial points
			this->set_ylm_conjg_data(iL, iM, thetas, phis, ylmData);

			SphericalPart integrand(f, c[0], c[1], c[2], ylmData);
			for (int iR = 0 ; iR < nR; ++iR )
			{
				double r = rgrid.get_radius(iR);
				integrand.set_radius(r);
				coefficientData[iR+nR*(Auxillary::memlayout::angular_momentum_layout(iL,iM))]
								= integrator.integrate(integrand, lebedev_rule);
			}
		}
	}

	this->initialize(lmax, std::move(coefficientData), std::move(rgrid));
}

inline Auxillary::alignedvector::ZV::iterator
SphericalHarmonicExpansion::begin()
{
	return data_.begin();
}

inline Auxillary::alignedvector::ZV::iterator
SphericalHarmonicExpansion::end()
{
	return data_.end();
}

inline Auxillary::alignedvector::ZV::const_iterator
SphericalHarmonicExpansion::begin() const
{
	return data_.begin();
}

inline Auxillary::alignedvector::ZV::const_iterator
SphericalHarmonicExpansion::end() const
{
	return data_.end();
}

} /* namespace AtomicSite */
} /* namespace elephon */
