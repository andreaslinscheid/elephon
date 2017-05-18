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
#include "IOMethods/ReadVASPSymmetries.h"

BOOST_AUTO_TEST_CASE( Build_Al_symmetries )
{
	boost::filesystem::path p(__FILE__);
	boost::filesystem::path dir = p.parent_path();
	std::string Al_test_outcar = std::string(dir.c_str())+"/../IOMethods/OUTCAR_Al_test.dat";

	elephon::IOMethods::ReadVASPSymmetries symReader;
	symReader.read_file(Al_test_outcar);

	elephon::LatticeStructure::Symmetry sym;
	sym.initialize( 1e-6, symReader.get_symmetries(), symReader.get_fractionTranslations() );

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
