/*	This file test_RegularGrid.cpp is part of elephon.
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
 *  Created on: May 23, 2017
 *      Author: A. Linscheid
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Input_test
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "IOMethods/VASPInterface.h"
#include "LatticeStructure/RegularGrid.h"
#include "fixtures/MockStartup.h"
#include <iostream>

BOOST_AUTO_TEST_CASE( Al_fcc_UnitCell )
{
	test::fixtures::MockStartup ms;
	auto testd = ms.get_data_for_testing_dir() / "Al" / "vasp" / "fcc_primitive";

	elephon::IOMethods::InputOptions noop;

	//introduce the VASP data loader
	std::shared_ptr<elephon::IOMethods::VASPInterface> loader =
			std::make_shared< elephon::IOMethods::VASPInterface >(noop);

	elephon::LatticeStructure::LatticeModule lattice;
	elephon::LatticeStructure::Symmetry sym;
	std::vector<elephon::LatticeStructure::Atom> atoms;
	elephon::LatticeStructure::RegularGrid kgrid;
	loader->read_cell_paramters( testd.string() ,1e-6,kgrid,lattice,atoms,sym);

	sym.set_reciprocal_space_sym();
	elephon::LatticeStructure::RegularGrid kgridDefaultIrredZone;
	kgridDefaultIrredZone.initialize(1e-6,
			kgrid.get_grid_dim(),
			kgrid.get_grid_shift(),
			sym,
			lattice);

	auto symIrredToRed = kgridDefaultIrredZone.get_maps_sym_irred_to_reducible();
	auto irredToRed = kgridDefaultIrredZone.get_maps_irreducible_to_reducible();
	auto symRedToIrred = kgridDefaultIrredZone.get_maps_sym_red_to_irreducible();
	auto redToIrred = kgridDefaultIrredZone.get_maps_red_to_irreducible();

	//Check that the mapping make sense and take one from the reducible to the irreducible k vector
	for ( int ikir = 0 ; ikir < kgridDefaultIrredZone.get_np_irred(); ++ikir)
	{
		auto kIrred = kgridDefaultIrredZone.get_vector_direct( irredToRed[ikir][0] );
		for ( int istar = 0 ; istar < int(symIrredToRed[ikir].size()); ++istar)
		{
			int reducibleIndex = irredToRed[ikir][istar];
			int irreducibleIndex = redToIrred[reducibleIndex];
			BOOST_REQUIRE( ikir == irreducibleIndex );

			int symIrRe = symIrredToRed[ikir][istar];
			int symReIr = symRedToIrred[reducibleIndex];
			BOOST_REQUIRE( (sym.get_index_inverse(symIrRe) == symReIr)
						&& (sym.get_index_inverse(symReIr) == symIrRe));
			auto kRed = kgridDefaultIrredZone.get_vector_direct( irredToRed[ikir][istar] );

			auto kIrredRot = kIrred;
			sym.apply(symIrRe,kIrredRot);
			for ( int i = 0 ; i < 3; ++i)
				BOOST_CHECK_CLOSE( kIrredRot[i] , kRed[i],  1e-6 );
		}
	}

	symIrredToRed = kgrid.get_maps_sym_irred_to_reducible();
	irredToRed = kgrid.get_maps_irreducible_to_reducible();
	symRedToIrred = kgrid.get_maps_sym_red_to_irreducible();
	redToIrred = kgrid.get_maps_red_to_irreducible();
	for ( int ikir = 0 ; ikir < kgrid.get_np_irred(); ++ikir)
	{
		auto kIrred = kgrid.get_vector_direct( irredToRed[ikir][0] );
		for ( int istar = 0 ; istar < int(symIrredToRed[ikir].size()); ++istar)
		{
			int reducibleIndex = irredToRed[ikir][istar];
			int irreducibleIndex = redToIrred[reducibleIndex];
			BOOST_REQUIRE( ikir == irreducibleIndex );

			int symIrRe = symIrredToRed[ikir][istar];
			int symReIr = symRedToIrred[reducibleIndex];
			BOOST_REQUIRE( (sym.get_index_inverse(symIrRe) == symReIr)
						&& (sym.get_index_inverse(symReIr) == symIrRe));
			auto kRed = kgrid.get_vector_direct( irredToRed[ikir][istar] );

			auto kIrredRot = kIrred;
			sym.apply(symIrRe,kIrredRot);
			for ( int i = 0 ; i < 3; ++i)
				BOOST_CHECK_CLOSE( kIrredRot[i] , kRed[i],  1e-6 );
		}
	}
}

BOOST_AUTO_TEST_CASE( FeSe_UnitCell )
{
	test::fixtures::MockStartup ms;
	auto testd = ms.get_data_for_testing_dir() / "FeSe"/"vasp"/"wfct"/"symmetric";

	elephon::IOMethods::InputOptions noop;

	//introduce the VASP data loader
	std::shared_ptr<elephon::IOMethods::VASPInterface> loader =
			std::make_shared< elephon::IOMethods::VASPInterface >(noop);

	elephon::LatticeStructure::LatticeModule lattice;
	elephon::LatticeStructure::Symmetry sym;
	std::vector<elephon::LatticeStructure::Atom> atoms;
	elephon::LatticeStructure::RegularGrid kgrid;
	loader->read_cell_paramters(testd.string(),1e-6,kgrid,lattice,atoms,sym);

	sym.set_reciprocal_space_sym();
	elephon::LatticeStructure::RegularGrid kgridDefaultIrredZone;
	kgridDefaultIrredZone.initialize(1e-6,
			kgrid.get_grid_dim(),
			kgrid.get_grid_shift(),
			sym,
			lattice);

	auto symIrredToRed = kgridDefaultIrredZone.get_maps_sym_irred_to_reducible();
	auto irredToRed = kgridDefaultIrredZone.get_maps_irreducible_to_reducible();
	auto symRedToIrred = kgridDefaultIrredZone.get_maps_sym_red_to_irreducible();
	auto redToIrred = kgridDefaultIrredZone.get_maps_red_to_irreducible();

	//Check that the mapping make sense and take one from the reducible to the irreducible k vector
	for ( int ikir = 0 ; ikir < kgridDefaultIrredZone.get_np_irred(); ++ikir)
	{
		auto kIrred = kgridDefaultIrredZone.get_vector_direct( irredToRed[ikir][0] );
		for ( int istar = 0 ; istar < int(symIrredToRed[ikir].size()); ++istar)
		{
			int reducibleIndex = irredToRed[ikir][istar];
			int irreducibleIndex = redToIrred[reducibleIndex];
			BOOST_REQUIRE( ikir == irreducibleIndex );

			int symIrRe = symIrredToRed[ikir][istar];
			int symReIr = symRedToIrred[reducibleIndex];
			BOOST_REQUIRE( (sym.get_index_inverse(symIrRe) == symReIr)
						&& (sym.get_index_inverse(symReIr) == symIrRe));
			auto kRed = kgridDefaultIrredZone.get_vector_direct( irredToRed[ikir][istar] );

			auto kIrredRot = kIrred;
			sym.apply(symIrRe,kIrredRot);
			for ( int i = 0 ; i < 3; ++i)
				BOOST_CHECK_CLOSE( kIrredRot[i] , kRed[i],  1e-6 );
		}
	}

	symIrredToRed = kgrid.get_maps_sym_irred_to_reducible();
	irredToRed = kgrid.get_maps_irreducible_to_reducible();
	symRedToIrred = kgrid.get_maps_sym_red_to_irreducible();
	redToIrred = kgrid.get_maps_red_to_irreducible();
	for ( int ikir = 0 ; ikir < kgrid.get_np_irred(); ++ikir)
	{
		auto kIrred = kgrid.get_vector_direct( irredToRed[ikir][0] );
		for ( int istar = 0 ; istar < int(symIrredToRed[ikir].size()); ++istar)
		{
			int reducibleIndex = irredToRed[ikir][istar];
			int irreducibleIndex = redToIrred[reducibleIndex];
			BOOST_REQUIRE( ikir == irreducibleIndex );

			int symIrRe = symIrredToRed[ikir][istar];
			int symReIr = symRedToIrred[reducibleIndex];
			BOOST_REQUIRE( (sym.get_index_inverse(symIrRe) == symReIr)
						&& (sym.get_index_inverse(symReIr) == symIrRe));
			auto kRed = kgrid.get_vector_direct( irredToRed[ikir][istar] );

			auto kIrredRot = kIrred;
			sym.apply(symIrRe,kIrredRot);
			for ( int i = 0 ; i < 3; ++i)
				BOOST_CHECK_CLOSE( kIrredRot[i] , kRed[i],  1e-6 );
		}
	}
}
