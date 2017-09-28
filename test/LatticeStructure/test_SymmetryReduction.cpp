/*	This file test_SymmetryReduction.cpp is part of elephon.
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
 *  Created on: Jun 18, 2017
 *      Author: A. Linscheid
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Input_test
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "LatticeStructure/SymmetryReduction.h"
#include "IOMethods/ReadVASPSymmetries.h"
#include "LatticeStructure/Atom.h"
#include "fixtures/MockStartup.h"
#include <vector>

BOOST_AUTO_TEST_CASE( irreducible_zones )
{
	using namespace elephon::LatticeStructure;
	typedef std::vector<double> V;
	std::vector<bool> F = {false,false,false};

	//Set up the set of reducible atoms in the conventional unit cell of Al
	std::vector<Atom> reducible(4, Atom(26.981, "Al",V({ 0.000000, 0.000000, 0.000000 }),F) );
	reducible[1] = Atom(26.981, "Al",V({ 0.000000, 0.500000, 0.500000 }),F);
	reducible[2] = Atom(26.981, "Al",V({ 0.500000, 0.000000, 0.500000 }),F);
	reducible[3] = Atom(26.981, "Al",V({ 0.500000, 0.500000, 0.000000 }),F);

	//load the symmetry
	elephon::IOMethods::ReadVASPSymmetries symReader;
	test::fixtures::MockStartup ms;
	auto testd = ms.get_data_for_testing_dir() / "Al" / "vasp" / "conventional";
	symReader.read_file( (testd / "OUTCAR").string() );
	elephon::LatticeStructure::LatticeModule lattice;
	lattice.initialize( std::vector<double>({1,0,0,0,1,0,0,0,1}) );
	elephon::LatticeStructure::Symmetry sym;
	sym.initialize( 1e-6, symReader.get_symmetries(), symReader.get_fractionTranslations(), lattice, true );

	//perform the reduction
	std::vector<Atom> irreducible;
	std::vector<int> redToIrred;
	std::vector<int> symRedToIrred;
	std::vector< std::vector<int> > irredToRed;
	std::vector< std::vector<int> > symIrredToRed;
	SymmetryReduction<Atom> symmReduce(sym, reducible,
			irreducible, redToIrred, symRedToIrred,
			irredToRed, symIrredToRed);

	for ( int i = 0 ; i < redToIrred.size(); ++i )
	{
		int irr = redToIrred[i];
		int isym = symRedToIrred[i];
		int isym_inv = sym.get_index_inverse(isym);
		std::map<int,int> star;
		for ( int istar = 0 ; istar < symIrredToRed[irr].size() ; ++istar)
			star.insert( std::make_pair(symIrredToRed[irr][istar],istar) );

		auto it = star.find(isym_inv);
		BOOST_REQUIRE( it != star.end() );

		int invstar = it->second;

		BOOST_REQUIRE( irredToRed[redToIrred[i]][invstar] == i );

		auto p = reducible[ i ].get_position();
		auto pr = irreducible[ irr ].get_position();

		sym.apply( symIrredToRed[irr][invstar], pr );

		BOOST_REQUIRE( p == pr );
	}
}
