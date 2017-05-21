/*	This file SymmetryReduction.hpp is part of elephon.
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
 *  Created on: May 21, 2017
 *      Author: A. Linscheid
 */

#include "SymmetryReduction.h"
#include <set>
#include <stdexcept>

namespace elephon
{
namespace LatticeStructure
{

template<class C>
SymmetryReduction<C>::SymmetryReduction()
{

}

template<class C>
SymmetryReduction<C>::SymmetryReduction(
		LatticeStructure::Symmetry const & sym,
		std::vector<C> const & reducible,
		std::vector<C> & irreducible,
		std::vector<int> & redToIrred,
		std::vector<int> & symRedToIrred,
		std::vector< std::vector<int> > & irredToRed,
		std::vector< std::vector<int> > & symIrredToRed )
{
	this->reduce(sym, reducible, irreducible, redToIrred, symRedToIrred, irredToRed, symIrredToRed);
}

template<class C>
void
SymmetryReduction<C>::reduce(
		LatticeStructure::Symmetry const & sym,
		std::vector<C> const & reducible,
		std::vector<C> & irreducible,
		std::vector<int> & redToIrred,
		std::vector<int> & symRedToIrred,
		std::vector< std::vector<int> > & irredToRed,
		std::vector< std::vector<int> > & symIrredToRed ) const
{
	//We apply the symmetry group and reduce the set to the non-symmetry-equivalent set
	//We also keep track of the index of the operations that maps back to the reducible set.
	//
	//	Algorithm:
	//		1) init redToIrred to values that do not appear, e.g. redToIrred.size()
	//		2) walk through the reducible array and insert the current index
	//			into the mapping of all _later_ indices that can be reached by a symmetry operation.
	//		3) skip indices that have already been mapped but record the connection
	//
	//	Result: Set of irreducible objects where we keep the objects that are
	//		_first_ in the list of reducible points.
	//
	//This requires a concept of 'equal', because we need to locate a transformed object.
	//For efficient look-up, in fact wee need and ordering by '<'.
	//This definition must be provided by the class C
	std::set<C> reducibleSet( reducible.begin(), reducible.end() );

	const int numRed = static_cast<int>(reducibleSet.size());
	const int numSym = sym.get_num_symmetries();
	redToIrred = std::vector<int> ( reducibleSet.size(), numRed );
	symRedToIrred = std::vector<int> ( reducibleSet.size(), numSym );

	//Note that we need irredIndex+symmetry index to go to the reducible set.
	//We don't know at this point how many irred displacements we have, so we init with numRed
	std::vector<std::vector<int>> indexIrreducibleToReducible(numRed,std::vector<int>(numSym)),
									symIrreducibleToReducible(numRed,std::vector<int>(numSym));

	std::vector<int> dimStarIrred(numRed,1);

	std::set<C> irreducibleSet;
	int idirr = 0;
	for (int id=0 ; id < numRed ; ++id)
	{
		//If the following is true, this reducible displacement was matched before -
		//do not add it to the set of irreducible ones but save it in the corresponding mappings
		if ( redToIrred[id] < numRed)
		{
			int ir = redToIrred[id];
			indexIrreducibleToReducible[ir][ dimStarIrred[ir] ] = id;
			symIrreducibleToReducible[ir][ dimStarIrred[ir] ] = symRedToIrred[id];
			dimStarIrred[ir]++;
			continue;
		}

		//insert this irreducible displacement into the set
		auto ret = irreducibleSet.insert( reducible[id] );
		if (not ret.second)
			throw std::logic_error("Trying to add a irreducible object that is already present");

		//record this displacement in the mappings
		int sid = sym.get_identity_index();
		indexIrreducibleToReducible[idirr][sid] = id;
		symIrreducibleToReducible[idirr][sid] = sid;

		//Rotate this displacement with all symmetry operators and find it in the reducible set.
		for (int isym=0;isym<numSym;isym++)
		{
			auto obj = reducible[id];
			obj.transform( sym.get_sym_op( isym ) );
			auto it = reducibleSet.find(obj);

			//The reducible grid is closed under its symmetry operations
			if ( it == reducibleSet.end() )
				throw std::logic_error("Reducible set of objects not "
						"closed under its symmetry operations");

			//check if this point was found before, but does not match the present irreducible index
			bool foundBefore = (redToIrred[id] != numRed);
			if ( foundBefore and (redToIrred[id] != idirr) )
				throw std::logic_error("Stars of objects are not distinct");

			//Add this irreducible index into all later points
			int indexInReducibleVector = std::distance(reducibleSet.begin(),it);
			if ( redToIrred[indexInReducibleVector] == numRed) {
				redToIrred[indexInReducibleVector] = idirr;
				symRedToIrred[indexInReducibleVector] = isym;
			}
		}
		idirr++;
	}

	irreducible = std::vector<C>(irreducibleSet.begin(),irreducibleSet.end());

	//Clean up - reduce the mappings to their actual size
	int numIrred = static_cast<int>(irreducible.size());
	irredToRed = std::vector< std::vector<int> >(numIrred);
	symIrredToRed = std::vector< std::vector<int> >(numIrred);
	for ( int ir = 0; ir < numIrred; ++ir)
	{
		irredToRed[ir] = std::vector<int>(dimStarIrred[ir]);
		symIrredToRed[ir] = std::vector<int>(dimStarIrred[ir]);
		for ( int isym = 0; isym < dimStarIrred[ir]; ++isym )
		{
			irredToRed[ir][isym] = indexIrreducibleToReducible[ir][isym];
			symIrredToRed[ir][isym] = symIrreducibleToReducible[ir][isym];
		}
	}
}

} /* namespace LatticeStructure */
} /* namespace elephon */
