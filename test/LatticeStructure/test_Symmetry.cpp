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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Input_test
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "LatticeStructure/Symmetry.h"
#include "LatticeStructure/LatticeModule.h"
#include "IOMethods/ReadVASPSymmetries.h"
#include "fixtures/MockStartup.h"

BOOST_AUTO_TEST_CASE( Build_Al_symmetries )
{
	test::fixtures::MockStartup ms;
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
	test::fixtures::MockStartup ms;
	auto testd = ms.get_data_for_testing_dir() / "Al" / "vasp" / "conventional";
	elephon::IOMethods::ReadVASPSymmetries symReader;
	symReader.read_file( (testd / "OUTCAR").string() );

	elephon::LatticeStructure::LatticeModule lattice;
	lattice.initialize( std::vector<double>({1,0,0,0,1,0,0,0,1}) );

	elephon::LatticeStructure::Symmetry sym;
	sym.initialize( 1e-6, symReader.get_symmetries(), symReader.get_fractionTranslations(), lattice, true );
}
