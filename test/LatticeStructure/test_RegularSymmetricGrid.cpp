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

#define BOOST_TEST_MODULE Input_test
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "IOMethods/VASPInterface.h"
#include "IOMethods/ResourceHandler.h"
#include <LatticeStructure/RegularSymmetricGrid.h>
#include "fixtures/MockStartup.h"
#include <iostream>

BOOST_AUTO_TEST_CASE( FeSe_vasp_k_points_file_reference )
{
	elephon::test::fixtures::MockStartup ms;
	auto testd = ms.get_data_for_testing_dir() / "FeSe" / "vasp" / "wfct" /"symmetric";

	std::string input = std::string()+
			"root_dir = "+testd.string()+"\n";
	elephon::IOMethods::InputOptions opts;
	ms.simulate_elephon_input(
			(testd / "infile").string(),
			input,
			opts);

	auto loader = std::make_shared<elephon::IOMethods::VASPInterface>(opts);

	elephon::LatticeStructure::RegularSymmetricGrid kgrid;
	loader->read_reciprocal_symmetric_grid(testd.string(), kgrid);

	// taken from the VASP symmetric output
	std::vector<double> known_irreducible_points{
		 0.125000, 0.125000, 0.250000 ,
		 0.375000, 0.125000, 0.250000 ,
		 0.375000, 0.375000, 0.250000 };
	std::vector<int> weights{ 8, 16, 8 };

	// taken from the VASP output with symmetry switched off
	std::vector<double> known_reducible_points{
		 0.125000, 0.125000, 0.250000  ,
		 0.375000, 0.125000, 0.250000 ,
		-0.375000, 0.125000, 0.250000 ,
		-0.125000, 0.125000, 0.250000 ,
		 0.125000, 0.375000, 0.250000 ,
		 0.375000, 0.375000, 0.250000 ,
		-0.375000, 0.375000, 0.250000 ,
		-0.125000, 0.375000, 0.250000 ,
		 0.125000,-0.375000, 0.250000 ,
		 0.375000,-0.375000, 0.250000 ,
		-0.375000,-0.375000, 0.250000 ,
		-0.125000,-0.375000, 0.250000 ,
		 0.125000,-0.125000, 0.250000 ,
		 0.375000,-0.125000, 0.250000 ,
		-0.375000,-0.125000, 0.250000 ,
		-0.125000,-0.125000, 0.250000 ,
		 0.125000, 0.125000,-0.250000 ,
		 0.375000, 0.125000,-0.250000 ,
		-0.375000, 0.125000,-0.250000 ,
		-0.125000, 0.125000,-0.250000 ,
		 0.125000, 0.375000,-0.250000 ,
		 0.375000, 0.375000,-0.250000 ,
		-0.375000, 0.375000,-0.250000 ,
		-0.125000, 0.375000,-0.250000 ,
		 0.125000,-0.375000,-0.250000 ,
		 0.375000,-0.375000,-0.250000 ,
		-0.375000,-0.375000,-0.250000 ,
		-0.125000,-0.375000,-0.250000 ,
		 0.125000,-0.125000,-0.250000 ,
		 0.375000,-0.125000,-0.250000 ,
		-0.375000,-0.125000,-0.250000 ,
		-0.125000,-0.125000,-0.250000  };

	BOOST_REQUIRE_EQUAL(kgrid.get_np_irred(), known_irreducible_points.size()/3 );

	BOOST_REQUIRE_EQUAL(kgrid.get_np_red(), known_reducible_points.size()/3 );

	for ( int ikirred = 0; ikirred  < kgrid.get_np_irred() ; ++ ikirred )
	{
		auto irtor = kgrid.get_maps_irreducible_to_reducible()[ikirred];
		BOOST_REQUIRE_EQUAL(irtor.size(), weights[ikirred]);
		auto symirtor = kgrid.get_maps_sym_irred_to_reducible()[ikirred];

		// confirm the irreducible k point coordinates
		int idSymOpIndex = kgrid.get_symmetry().get_identity_index();
		auto irredvec = kgrid.get_vector_direct(irtor[idSymOpIndex]);
		std::vector<double> ref_irredvec( &known_irreducible_points[ikirred*3],
										  &known_irreducible_points[ikirred*3]+3);
		for ( int i = 0; i < 3; ++i)
			BOOST_CHECK_SMALL(irredvec[i]-ref_irredvec[i], kgrid.get_grid_prec());

		// confirm the reducible k point connection
		for ( int istar = 0 ; istar < symirtor.size() ; ++istar)
		{
			int ikred = irtor[istar];
			auto redvec = kgrid.get_vector_direct(ikred);
			std::vector<double> ref_redvec( &known_reducible_points[ikred*3],
										    &known_reducible_points[ikred*3]+3);
			for ( int i = 0; i < 3; ++i)
				BOOST_CHECK_SMALL(redvec[i]-ref_redvec[i], kgrid.get_grid_prec());
		}
	}
}

BOOST_AUTO_TEST_CASE( Al_fcc_UnitCell )
{
	elephon::test::fixtures::MockStartup ms;
	auto testd = ms.get_data_for_testing_dir() / "Al" / "vasp" / "fcc_primitive";

	elephon::IOMethods::InputOptions noop;

	//introduce the VASP data loader
	std::shared_ptr<elephon::IOMethods::VASPInterface> loader =
			std::make_shared< elephon::IOMethods::VASPInterface >(noop);

	elephon::LatticeStructure::LatticeModule lattice;
	elephon::LatticeStructure::Symmetry sym;
	std::vector<elephon::LatticeStructure::Atom> atoms;
	elephon::LatticeStructure::RegularSymmetricGrid kgrid;
	loader->read_cell_paramters( testd.string() ,1e-6,kgrid,lattice,atoms,sym);

	sym.set_reciprocal_space_sym();
	elephon::LatticeStructure::RegularSymmetricGrid kgridDefaultIrredZone;
	kgridDefaultIrredZone.initialize(
			kgrid.get_grid_dim(),
			1e-6,
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
	elephon::test::fixtures::MockStartup ms;
	auto testd = ms.get_data_for_testing_dir() / "FeSe"/"vasp"/"wfct"/"symmetric";

	elephon::IOMethods::InputOptions noop;

	//introduce the VASP data loader
	std::shared_ptr<elephon::IOMethods::VASPInterface> loader =
			std::make_shared< elephon::IOMethods::VASPInterface >(noop);

	elephon::LatticeStructure::LatticeModule lattice;
	elephon::LatticeStructure::Symmetry sym;
	std::vector<elephon::LatticeStructure::Atom> atoms;
	elephon::LatticeStructure::RegularSymmetricGrid kgrid;
	loader->read_cell_paramters(testd.string(),1e-6,kgrid,lattice,atoms,sym);

	sym.set_reciprocal_space_sym();
	elephon::LatticeStructure::RegularSymmetricGrid kgridDefaultIrredZone;
	kgridDefaultIrredZone.initialize(
			kgrid.get_grid_dim(),
			1e-6,
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

BOOST_AUTO_TEST_CASE( MgB2_k_points )
{
	elephon::test::fixtures::MockStartup ms;
	auto testd = ms.get_data_for_testing_dir() / "MgB2" / "vasp" / "ldos";

	std::string input = std::string()+
			"root_dir = "+testd.string()+"\n";
	elephon::IOMethods::InputOptions opts;
	ms.simulate_elephon_input(
			(testd / "infile").string(),
			input,
			opts);

	auto loader = std::make_shared<elephon::IOMethods::VASPInterface>(opts);

	elephon::LatticeStructure::RegularSymmetricGrid kgrid;
	elephon::LatticeStructure::LatticeModule lattice;
	std::vector<elephon::LatticeStructure::Atom> atoms;
	elephon::LatticeStructure::Symmetry sym;
	loader->read_cell_paramters(
			loader->get_optns().get_root_dir(),
			loader->get_optns().get_gPrec(),
			kgrid,
			lattice,
			atoms,
			sym);

	//The methodology is that we reproduce the irreducible zone of VASP exactly
	//thus, we take from the output file, the coordinates and the weights:
	std::vector<double> ref_irred_kpts_vec{
		0.000000, 0.000000, 0.000000,
		0.333333, 0.000000, 0.000000,
		0.333333, 0.333333, 0.000000,
		0.000000, 0.000000, 0.500000,
		0.333333, 0.000000, 0.500000,
		0.333333, 0.333333, 0.500000
	};
	std::vector<int> ref_k_vect_weights{1, 6, 2, 1, 6, 2};

	int nptotal = 0;
	for ( auto w : ref_k_vect_weights )
		nptotal += w;

	BOOST_CHECK_EQUAL(nptotal, 3*3*2);

	BOOST_CHECK_EQUAL(kgrid.get_np_red(), nptotal);

	BOOST_CHECK_EQUAL(kgrid.get_np_irred(), ref_k_vect_weights.size());

	for ( int ikirred = 0; ikirred  < kgrid.get_np_irred(); ++ ikirred )
	{
		auto irtor = kgrid.get_maps_irreducible_to_reducible()[ikirred];
		auto symirtor = kgrid.get_maps_sym_irred_to_reducible()[ikirred];
		//confirm that the multiplicity of each k point is correct
		BOOST_CHECK_EQUAL(irtor.size(), ref_k_vect_weights[ikirred]);

		//The first element in the star is the irreducible vector. Confirm
		//that this agrees with the zone used in VASP while taking into account the different conventions
		auto irredvec = kgrid.get_vector_direct(irtor[0]);
		std::vector<double> ref_irredvec(&ref_irred_kpts_vec[ikirred*3], &ref_irred_kpts_vec[ikirred*3]+3);
		for ( int i = 0; i < 3; ++i)
		{
			ref_irredvec[i] -= std::floor(ref_irredvec[i]+0.5);
			BOOST_CHECK_SMALL(irredvec[i]-ref_irredvec[i], kgrid.get_grid_prec());
		}

		for ( int is = 1; is < symirtor.size(); ++is )
		{
			int ikred = irtor[is];
			auto redvec_by_index = kgrid.get_vector_direct(ikred);
			int isym = symirtor[is];
			auto redvec = irredvec;
			kgrid.get_symmetry().apply(isym,redvec,true);
			for ( int i = 0; i < 3; ++i)
				BOOST_CHECK_SMALL(redvec[i]-redvec_by_index[i], kgrid.get_grid_prec());
		}
	}
}

