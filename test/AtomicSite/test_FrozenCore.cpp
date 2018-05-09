/*	This file test_FrozenCore.cpp is part of elephon.
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
 *  Created on: May 8, 2018
 *      Author: A. Linscheid
 */
#include <boost/test/unit_test.hpp>
#include "AtomicSite/FrozenCore.h"
#include "AtomicSite/RadialGrid.h"
#include "fixtures/DataLoader.h"
#include <vector>
#include <memory>

BOOST_AUTO_TEST_SUITE( FrozenCore )

/**
 * Create a test object that has a core charge of -Z and a homogeneous electronic charge
 * centered around the center with radius r. The charge integrates to e.
 *
 * @param[in] Z		Point charge of the core
 * @param[in] e		Integral of the electornic charge
 * @param[in] r		radial extend of the homogenous electronic charge
 * @return a shared pointer with the initialized object
 */
std::shared_ptr<elephon::AtomicSite::FrozenCore>
create_test_obj_homog_sphere_charge(double Z, double e, double r)
{
	auto fc = std::make_shared<elephon::AtomicSite::FrozenCore>();
	elephon::test::fixtures::DataLoader dl;
	elephon::AtomicSite::RadialGrid rgrid = dl.radial_sample_rgrid();
	std::vector<double> chargeTimesR2(rgrid.get_num_R());
	if (rgrid.get_range_of_definition() < r )
		throw std::logic_error("Choose an r smaller than the radius of the sample grid for testing.");
	double totalChargeVolume = 0.0;
	for (int ir = 0 ; ir < rgrid.get_num_R(); ++ir)
	{
		if (rgrid.get_radius(ir)<=r)
			totalChargeVolume = std::max(totalChargeVolume,4*M_PI*std::pow(rgrid.get_radius(ir),3)/3.0);
		chargeTimesR2[ir] = rgrid.get_radius(ir)<=r ? rgrid.get_radius(ir)*rgrid.get_radius(ir) : 0.0;
	}
	for (auto &c : chargeTimesR2)
		c *= e/totalChargeVolume;

	fc->initialize(Z,std::move(chargeTimesR2), rgrid);
	return fc;
}

BOOST_AUTO_TEST_CASE( test_charge )
{
	const double eleChg = 1.0;
	auto fc = create_test_obj_homog_sphere_charge(2.0,eleChg,0.5);
	BOOST_CHECK_SMALL(fc->total_electronic_charge()-eleChg, 5e-2);
}

BOOST_AUTO_TEST_CASE( test_potential )
{
	double eleChg = 1.0;
	double Z = 1.0;
	const double radiusCharge = 0.5;
	auto fc = create_test_obj_homog_sphere_charge(Z,eleChg,radiusCharge);

	// since the simpson rule integration is not well adopted to the step function
	// cutoff of this model some inaccuracy in the charge is unavoidable.
	eleChg = fc->total_electronic_charge();

	elephon::AtomicSite::RadialGrid const & rgrid = fc->get_radial_grid();
	const int nR = rgrid.get_num_R();

	// test the -Z/r potential first.
	std::vector<double> potential(nR, 0.0);
	fc->add_core_potential(potential.begin(), potential.end());
	double diffAnalytic = 0.0;
	for (int ir = 0 ; ir < nR; ++ir)
	{
		double analytic = -elephon::Auxillary::units::HARTREE_TO_EV*Z/rgrid.get_radius(ir)*std::sqrt(4*M_PI);
		diffAnalytic += std::abs( analytic - potential[ir]);
	}
	BOOST_CHECK_SMALL(diffAnalytic,1e-5);

	diffAnalytic = 0.0;
	std::fill(potential.begin(), potential.end(), 0.0);
	fc->add_core_hartree_potential(potential.begin(), potential.end());
	for (int ir = 0 ; ir < nR; ++ir)
	{
		double forzenElectronChargeLower = elephon::Auxillary::units::HARTREE_TO_EV*eleChg*std::sqrt(4*M_PI)
										*(3.0*std::pow(radiusCharge,2)-std::pow(rgrid.get_radius(ir),2))
										/2.0/std::pow(radiusCharge,3);
		double forzenElectronChargeGreater = elephon::Auxillary::units::HARTREE_TO_EV*eleChg/rgrid.get_radius(ir)*std::sqrt(4*M_PI);
		double analytic = rgrid.get_radius(ir) <= radiusCharge ? forzenElectronChargeLower:	forzenElectronChargeGreater;
		double dr = ir < nR-1 ? rgrid.get_radius(ir+1) -rgrid.get_radius(ir) : rgrid.get_radius(ir)-rgrid.get_radius(ir-1);
		diffAnalytic += std::abs((analytic-potential[ir])/(analytic+potential[ir]))*dr;
	}
	BOOST_CHECK_SMALL(diffAnalytic,5e-2);
}

template<class analyticRef>
void check_analytic_ref_on_grid(
		std::shared_ptr<elephon::AtomicSite::FrozenCore> fc,
		elephon::AtomicSite::RadialGrid const & rgrid,
		analyticRef method,
		double tolerance)
{
	const int nPtsCheckGrid = 30;

	elephon::test::fixtures::DataLoader::RegularGridAtom regGridAtom(nPtsCheckGrid, rgrid);

	elephon::AtomicSite::SphericalHarmonicExpansion she;
	she.initialize(1,elephon::Auxillary::alignedvector::ZV(4*rgrid.get_num_R()), rgrid);
	fc->add_core_hartree_analytic_displacement({0.0, 0.0, 1.0}, she.begin(), she.end());
	elephon::Auxillary::alignedvector::ZV dataZ, dataX, dataY;
	she.interpolate(regGridAtom.get_atom_grid(), dataZ);

	std::fill(she.begin(), she.end(), std::complex<double>(0.0));
	fc->add_core_hartree_analytic_displacement({0.0, 1.0, 0.0}, she.begin(), she.end());
	she.interpolate(regGridAtom.get_atom_grid(), dataY);

	std::fill(she.begin(), she.end(), std::complex<double>(0.0));
	fc->add_core_hartree_analytic_displacement({1.0, 0.0, 0.0}, she.begin(), she.end());
	she.interpolate(regGridAtom.get_atom_grid(), dataX);

	int c = 0;
	double diffReal = 0.0, diffImag = 0.0;
	double x,y,z;
	for (int ix = 0 ; ix < nPtsCheckGrid; ++ix)
		for (int iy = 0 ; iy < nPtsCheckGrid; ++iy)
			for (int iz = 0 ; iz < nPtsCheckGrid; ++iz)
				if ( regGridAtom.check_coord_in_atomic_sphere(ix,iy,iz,x,y,z) )
				{
					double referenceValue = method(x,y,z,z);
					diffReal += std::abs(std::real(dataZ[c])-referenceValue);
					diffImag += std::abs(std::imag(dataZ[c]));

					referenceValue = method(x,y,z,y);
					diffReal += std::abs(std::real(dataY[c])-referenceValue);
					diffImag += std::abs(std::imag(dataY[c]));

					referenceValue = method(x,y,z,x);
					diffReal += std::abs(std::real(dataX[c])-referenceValue);
					diffImag += std::abs(std::imag(dataX[c]));
					++c;
				}

	diffReal /= regGridAtom.get_atom_grid().size();
	diffImag /= regGridAtom.get_atom_grid().size();
	BOOST_CHECK_SMALL(diffReal, tolerance);
	BOOST_CHECK_SMALL(diffImag, tolerance);
}

BOOST_AUTO_TEST_CASE( test_potential_displacement )
{
	double eleChg = 0.0; // first do just core point charge displacement
	double Z = 1.0;
	const double radiusCharge = 0.5;
	auto fc = create_test_obj_homog_sphere_charge(Z,eleChg,radiusCharge);

	// test for the classical dipole
	auto dipoleDir = [&] (double x, double y, double z, double dir) {
		double r = std::sqrt(std::pow(x,2)+std::pow(y,2)+std::pow(z,2));
		return - elephon::Auxillary::units::HARTREE_TO_EV*Z /(r*r*r) * dir;
	};

	check_analytic_ref_on_grid(fc, fc->get_radial_grid(), dipoleDir, 1e-3); // analytic, thus quite accurate

	// now check the frozen core charge only
	eleChg = 1.0;
	Z = 0.0;
	fc = create_test_obj_homog_sphere_charge(Z,eleChg,radiusCharge);
	eleChg = fc->total_electronic_charge();
	double density = eleChg / (4*M_PI/3.0*std::pow(radiusCharge,3));

	// test for the classical dipole
	auto dipoleFrozenEleChg = [&] (double x, double y, double z, double dir) {
		double r = std::sqrt(std::pow(x,2)+std::pow(y,2)+std::pow(z,2));
		if ( r > radiusCharge)
			return elephon::Auxillary::units::HARTREE_TO_EV*eleChg /(r*r*r) * dir;
		return elephon::Auxillary::units::HARTREE_TO_EV * dir / 3.0 * (4*M_PI*density);
	};
	check_analytic_ref_on_grid(fc, fc->get_radial_grid(), dipoleFrozenEleChg, 0.5); // numeric, hard cutoff, large values => not very accurate
}

BOOST_AUTO_TEST_SUITE_END()
