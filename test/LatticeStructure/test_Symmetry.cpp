/*	This file test_Symmetry.cpp is part of elephon.
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
 *  Created on: May 16, 2017
 *      Author: A. Linscheid
 */
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "LatticeStructure/Symmetry.h"
#include "LatticeStructure/LatticeModule.h"
#include "IOMethods/VASPInterface.h"
#include "IOMethods/ReadVASPSymmetries.h"
#include "fixtures/MockStartup.h"

BOOST_AUTO_TEST_SUITE( Symmetry )

BOOST_AUTO_TEST_CASE( Si_vasp )
{
	elephon::test::fixtures::MockStartup ms;
	auto testd = ms.get_data_for_testing_dir() / "Si" / "vasp" / "conventional";
	std::string input = std::string()+
			"root_dir = "+testd.string()+"\n";
	elephon::IOMethods::InputOptions opts;
	ms.simulate_elephon_input(
			(testd / "infile").string(),
			input,
			opts);

	auto loader = std::make_shared<elephon::IOMethods::VASPInterface>(opts);
	elephon::LatticeStructure::LatticeModule lat;
	loader->read_lattice_structure(opts.get_root_dir(), lat);

	elephon::LatticeStructure::Symmetry sym;
	loader->read_symmetry( opts.get_root_dir(), 1e-4, lat, sym);

	BOOST_CHECK( sym.get_num_symmetries() == 48 );
	BOOST_CHECK( sym.has_inversion() == false);
}

BOOST_AUTO_TEST_CASE( test_MgB2 )
{
	elephon::test::fixtures::MockStartup ms;
	auto testd = ms.get_data_for_testing_dir() / "MgB2" / "vasp" / "ldos";

	elephon::IOMethods::InputOptions opts;
	ms.simulate_elephon_input(
			(testd/"infile").string(),
			std::string("root_dir=")+testd.string()+"\n",
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

	BOOST_REQUIRE( sym.get_num_symmetries() == 24 );
}

BOOST_AUTO_TEST_CASE( Build_Al_symmetries )
{
	elephon::test::fixtures::MockStartup ms;
	auto testd = ms.get_data_for_testing_dir() / "Al" / "vasp" / "conventional";
	elephon::IOMethods::ReadVASPSymmetries symReader;
	symReader.read_file( (testd / "OUTCAR").string() );

	elephon::LatticeStructure::LatticeModule lattice;
	lattice.initialize( std::vector<double>({1,0,0,0,1,0,0,0,1}) );

	elephon::LatticeStructure::Symmetry sym;
	sym.initialize( 1e-6, symReader.get_symmetries(), symReader.get_fractionTranslations(), lattice, true );

	BOOST_REQUIRE( sym.get_num_symmetries() == 48 );

	int refIndex  = 22;
	std::vector<int> ref_sym_22 = {0,0,1,0,-1,0,1,0,0};
	auto sop = sym.get_sym_op(refIndex);
	for ( int i = 0; i < 3 ; ++i)
	{
		for ( int j = 0; j < 3 ; ++j)
			BOOST_REQUIRE(sop.ptgroup[i*3+j] == ref_sym_22[i*3+j]);
		BOOST_REQUIRE( std::fabs(sop.fracTrans[i]) < 1e-16);
	}

	//We check the multiplication of irot 4 with 17 which matches irot 41
	//	/0	1	0\	/0	1	0\ 		0	0	1
	//	|0	0	1|*	|0	0	1| =  	-1	0	0
	//	\1	0	0/	\-1	0	0/ 		0	1	0
	std::vector<int> ref_sym_4 = 	{0,1,0, 0,0,1,  1,0,0};
	std::vector<int> ref_sym_17 = 	{0,1,0, 0,0,1, -1,0,0};
	std::vector<int> ref_sym_41 = 	{0,0,1, -1,0,0, 0,1,0};
	auto sop4 = sym.get_sym_op(4);
	auto sop17 = sym.get_sym_op(17);
	auto sop41 = sym.get_sym_op(41);
	for ( int i = 0; i < 3 ; ++i)
		for ( int j = 0; j < 3 ; ++j)
		{
			BOOST_REQUIRE(sop4.ptgroup[i*3+j] == ref_sym_4[i*3+j]);
			BOOST_REQUIRE(sop17.ptgroup[i*3+j] == ref_sym_17[i*3+j]);
			BOOST_REQUIRE(sop41.ptgroup[i*3+j] == ref_sym_41[i*3+j]);
		}
	BOOST_REQUIRE(sym.get_group_product(4,17) == 41 );

	//Confirm that irot 26 has the inverse rot 40
	//		 0   1   0      0	 0	-1
	//26 =	 0   0  -1	40=	1	 0	 0
	//		-1   0   0 		0	-1	 0
	BOOST_REQUIRE(sym.get_index_inverse(26) == 40 );


	//Now we discard everything but inversion and identity
	std::vector<int> dopIndices;
	for ( int i = 2; i < sym.get_num_symmetries() ; ++i)
		dopIndices.push_back(i);
	sym.symmetry_reduction( dopIndices );
	BOOST_REQUIRE( sym.get_num_symmetries() == 2 );
	std::vector<int> ref_sym_2 = {-1,0,0,0,-1,0,0,0,-1};
	auto sopInv = sym.get_sym_op(1);
	for ( int i = 0; i < 3 ; ++i)
	{
		for ( int j = 0; j < 3 ; ++j)
			BOOST_REQUIRE(sopInv.ptgroup[i*3+j] == ref_sym_2[i*3+j]);
		BOOST_REQUIRE( std::fabs(sopInv.fracTrans[i]) < 1e-16);
	}
}

BOOST_AUTO_TEST_CASE( Symemtry_reduction )
{
	elephon::test::fixtures::MockStartup ms;
	auto testd = ms.get_data_for_testing_dir() / "Al" / "vasp" / "conventional";
	elephon::IOMethods::ReadVASPSymmetries symReader;
	symReader.read_file( (testd / "OUTCAR").string() );

	elephon::LatticeStructure::LatticeModule lattice;
	lattice.initialize( std::vector<double>({1,0,0,0,1,0,0,0,1}) );

	elephon::LatticeStructure::Symmetry sym;
	sym.initialize( 1e-6, symReader.get_symmetries(), symReader.get_fractionTranslations(), lattice, true );
}

BOOST_AUTO_TEST_CASE( MgB2_vasp )
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
	elephon::LatticeStructure::LatticeModule lat;
	std::vector<elephon::LatticeStructure::Atom> atoms;
	elephon::LatticeStructure::Symmetry sym;
	loader->read_cell_paramters(
			testd.string(),
			loader->get_optns().get_gPrec(),
			kgrid,
			lat,
			atoms,
			sym);

	//Check that all k vectors in the star have the same length
	auto itor = kgrid.get_maps_irreducible_to_reducible();
	for ( int ik = 0 ; ik < kgrid.get_np_irred(); ++ik )
	{
		auto kvecirrd = kgrid.get_vector_direct( itor[ik][0] );
		auto kvecIcart = kvecirrd;
		kgrid.get_lattice().reci_direct_to_cartesian(kvecIcart);
		double ni2 = std::sqrt(kvecIcart[0]*kvecIcart[0]+kvecIcart[1]*kvecIcart[1]+kvecIcart[2]*kvecIcart[2]);
		for ( int istar = 0 ; istar < itor[ik].size(); ++istar )
		{
			auto kvecread = kgrid.get_vector_direct( itor[ik][istar] );
			auto kvecRcart = kvecread;
			kgrid.get_lattice().reci_direct_to_cartesian(kvecRcart);
			double nr2 = std::sqrt(kvecRcart[0]*kvecRcart[0]+kvecRcart[1]*kvecRcart[1]+kvecRcart[2]*kvecRcart[2]);
			BOOST_CHECK_SMALL( ni2-nr2, 0.00001 );
		}
	}
}
BOOST_AUTO_TEST_SUITE_END()
