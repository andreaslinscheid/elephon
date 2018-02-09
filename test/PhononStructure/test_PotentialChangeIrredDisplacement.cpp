/*	This file test_PotentialChangeIrredDisplacement.cpp is part of elephon.
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
 *  Created on: Feb 3, 2018
 *      Author: A. Linscheid
 */
#include <boost/test/unit_test.hpp>
#include "PhononStructure/PotentialChangeIrredDisplacement.h"
#include "LatticeStructure/RegularBareGrid.h"
#include "AtomicSite/AtomSiteData.h"
#include "LatticeStructure/PrimitiveToSupercellConnection.h"
#include "LatticeStructure/UnitCell.h"
#include "LatticeStructure/Atom.h"
#include <memory>
#include <numeric>

BOOST_AUTO_TEST_SUITE( PotentialChangeIrredDisplacement )

BOOST_AUTO_TEST_CASE( test_explicit_example )
{
	using namespace elephon;
	// The potential change is the difference between two potential objects
	// on two different grids, the regular and the atomic site grid.
	// Here, we test the class
	// 1) by setting the data to be the same in the primitive cell part
	// 	so that the difference must be zero
	// 2) and to a known difference so that the expected integral over that area is known analytically.
	//  We also set the data in the radial part to zero in one instance in a given channel so that the data
	//	in that channel is purely from the other object.
	const int lMax = 5;
	const double Radius = 0.5;
	const int nRad = 10;

	auto ucGrid_ptr = std::make_shared<LatticeStructure::RegularBareGrid>(std::move(LatticeStructure::RegularBareGrid({10, 10, 10})));
	auto scGrid_ptr = std::make_shared<LatticeStructure::RegularBareGrid>(std::move(LatticeStructure::RegularBareGrid({30, 30, 30})));
	auto primCell_ptr = std::make_shared<LatticeStructure::UnitCell>();
	primCell_ptr->initialize(
			std::vector<LatticeStructure::Atom>(1, LatticeStructure::Atom(1.0, "A", {0.0, 0.0, 0.0}, {false, false, false})),
			LatticeStructure::LatticeModule(),
			LatticeStructure::Symmetry());
	auto superCell_ptr = std::make_shared<LatticeStructure::UnitCell>(std::move(primCell_ptr->build_supercell(3,3,3)));

	auto primScCon = std::make_shared<LatticeStructure::PrimitiveToSupercellConnection>();
	primScCon->initialize(primCell_ptr, superCell_ptr);

	// here we set the ground state "potential" on the regular grid to Pi
	std::vector<double> regularGridGroundStatePotential(ucGrid_ptr->get_num_points(), M_PI);

	// here we set the displaced "potential" on the regular grid to Pi
	std::vector<double> regularGridDisplacedPotential(scGrid_ptr->get_num_points(), M_PI);

	std::vector<double> rpoints(nRad);
	for (int ir = 0 ; ir < nRad; ++ir)
		rpoints[ir] = std::pow(static_cast<double>(ir+1)/static_cast<double>(nRad+1), 2);
	AtomicSite::RadialGrid rGrid;
	rGrid.initialize(primCell_ptr->get_atoms_list().front().get_position(), Radius, std::move(rpoints));

	// set the rotationally invariant part to 1;
	AtomicSite::SphericalHarmonicExpansion sexp;
	sexp.initialize(lMax, Auxillary::alignedvector::ZV((lMax+1)*(lMax+1)*nRad, 0.0), std::move(rGrid));
	for (int ir = 0; ir < sexp.get_radial_grid().get_num_R(); ++ir)
		sexp(ir, 0, 0) = std::complex<double>(1.0);

	AtomicSite::AtomSiteData adata;
	adata.initialize(
			primCell_ptr->get_atoms_list().front(),
			sexp);
	auto radialGroundStatePotential = std::vector<std::shared_ptr<const AtomicSite::AtomSiteData>>(1);
	for (auto &rgp_ptr : radialGroundStatePotential)
		rgp_ptr = std::make_shared<AtomicSite::AtomSiteData>(adata);
	auto radialDisplacedPotential = std::vector<std::shared_ptr<const AtomicSite::AtomSiteData>>(primScCon->supercell_volume_factor());

	// create the displaced atom and set the first element
	auto displacedAtom = primCell_ptr->get_atoms_list().front();
	auto displ_ptr = std::make_shared<LatticeStructure::AtomDisplacement>(displacedAtom, 0.001, std::vector<double>({1.0, 0.0, 0.0}));
	displacedAtom.apply_displacement(*displ_ptr);
	adata.initialize(displacedAtom, sexp);
	radialDisplacedPotential[0] = std::make_shared<AtomicSite::AtomSiteData>(std::move(adata));
	// set the remainder non-displaced atoms
	for (int ia = 1 ; ia < primScCon->supercell_volume_factor(); ++ia)
	{
		AtomicSite::AtomSiteData adata;
		adata.initialize(superCell_ptr->get_atoms_list()[ia], sexp);
		radialDisplacedPotential[ia] = std::make_shared<AtomicSite::AtomSiteData>(std::move(adata));
	}

	PhononStructure::PotentialChangeIrredDisplacement chng;
	chng.initialize(
			displ_ptr,
			regularGridGroundStatePotential,
			radialGroundStatePotential,
			regularGridDisplacedPotential,
			radialDisplacedPotential,
			ucGrid_ptr,
			scGrid_ptr,
			primScCon);

	// Check that the entire data in the object is zero numerically
	double sum = std::accumulate(chng.begin_regular_data(), chng.end_regular_data(), 0.0,
			[](double const & a, double const & b){return std::abs(a)+std::abs(b);});
	BOOST_CHECK_SMALL(sum, 1e-8);

	for (int iA = 0; iA < superCell_ptr->get_atoms_list().size(); ++iA)
	{
		std::complex<double> sum_c = std::accumulate(chng.begin_radial_data(iA), chng.end_radial_data(iA), std::complex<double>(0.0),
				[](std::complex<double> const & a, std::complex<double> const & b){return std::abs(a)+std::abs(b);});
		BOOST_CHECK_SMALL(std::abs(sum_c), 1e-8);
	}

	// now point 2)
	std::fill(regularGridDisplacedPotential.begin(), regularGridDisplacedPotential.end(), 2.0*M_PI);

	// create the displaced atom and set the first element. Here we don't actually displace it, for the following reason:
	// We want to test the data generation from the finite difference with a known value of the spherical harmonic coefficient.
	// We thus simply assume the data fitting for a shifted position is working as expected (see tests of AtomicSite::SphericalHarmonicExpansion).
	const int lCheck = 2;
	const int mCheck = -1;
	const std::complex<double> zeroLM = std::complex<double>(2.0);
	const std::complex<double> finiteLM = std::complex<double>(3.0);
	for (int ir = 0; ir < sexp.get_radial_grid().get_num_R(); ++ir)
	{
		sexp(ir, 0, 0) = zeroLM;
		sexp(ir, mCheck, lCheck) = finiteLM;
	}
	adata.initialize(primCell_ptr->get_atoms_list().front(), sexp);
	radialDisplacedPotential[0] = std::make_shared<AtomicSite::AtomSiteData>(adata);
	const int lastElem = static_cast<int>(radialDisplacedPotential.size()) - 1;
	adata.initialize(radialDisplacedPotential[lastElem]->get_atom(), sexp);
	radialDisplacedPotential[lastElem] = std::make_shared<AtomicSite::AtomSiteData>(adata);

	chng.initialize(
			displ_ptr,
			regularGridGroundStatePotential,
			radialGroundStatePotential,
			regularGridDisplacedPotential,
			radialDisplacedPotential,
			ucGrid_ptr,
			scGrid_ptr,
			primScCon);

	sum = std::accumulate(chng.begin_regular_data(), chng.end_regular_data(), 0.0,
			[](double const & a, double const & b){return std::abs(a)+std::abs(b);});
	BOOST_CHECK_CLOSE(sum, M_PI*regularGridDisplacedPotential.size(), 1e-8);

	// check first atom
	for (int il = 0; il <= lMax; ++il)
		for (int im = -il; im <= il; ++im)
		{
			auto itB = chng.begin_radial_data(0)+Auxillary::memlayout::angular_momentum_layout(il,im)*nRad;
			auto itE = itB + nRad;
			std::complex<double> sum_c = std::accumulate(itB, itE, std::complex<double>(0.0),
					[](std::complex<double> const & a, std::complex<double> const & b){return a+b;});
			if ((il == 0) && (im == 0))
				BOOST_CHECK_SMALL(std::abs(sum_c - (std::complex<double>(1.0)-zeroLM)*static_cast<double>(nRad)), 1e-8);
			else if ((il == lCheck) && (im == mCheck))
				BOOST_CHECK_SMALL(std::abs(sum_c - (-finiteLM)*static_cast<double>(nRad)), 1e-8);
			else
				BOOST_CHECK_SMALL(std::abs(sum_c), 1e-8);
		}

	// check last atom
	for (int il = 0; il <= lMax; ++il)
		for (int im = -il; im <= il; ++im)
		{
			auto itB = chng.begin_radial_data(lastElem)+Auxillary::memlayout::angular_momentum_layout(il,im)*nRad;
			auto itE = itB + nRad;
			std::complex<double> sum_c = std::accumulate(itB, itE, std::complex<double>(0.0),
					[](std::complex<double> const & a, std::complex<double> const & b){return a+b;});
			if ((il == 0) && (im == 0))
				BOOST_CHECK_SMALL(std::abs(sum_c - (std::complex<double>(1.0)-zeroLM)*static_cast<double>(nRad)), 1e-8);
			else if ((il == lCheck) && (im == mCheck))
				BOOST_CHECK_SMALL(std::abs(sum_c - (-finiteLM)*static_cast<double>(nRad)), 1e-8);
			else
				BOOST_CHECK_SMALL(std::abs(sum_c), 1e-8);
		}

	// check all the other atoms
	for (int iA = 1; iA < superCell_ptr->get_atoms_list().size(); ++iA)
	{
		if ( iA == lastElem )
			continue;
		std::complex<double> sum_c = std::accumulate(chng.begin_radial_data(iA), chng.end_radial_data(iA), std::complex<double>(0.0),
				[](std::complex<double> const & a, std::complex<double> const & b){return std::abs(a)+std::abs(b);});
		BOOST_CHECK_SMALL(std::abs(sum_c), 1e-8);
	}
}

BOOST_AUTO_TEST_SUITE_END()
