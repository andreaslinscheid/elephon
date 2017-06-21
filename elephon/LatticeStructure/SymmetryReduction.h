/*	This file SymmetryReduction.h is part of elephon.
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

#ifndef ELEPHON_LATTICESTRUCTURE_SYMMETRYREDUCTION_H_
#define ELEPHON_LATTICESTRUCTURE_SYMMETRYREDUCTION_H_

#include "LatticeStructure/Symmetry.h"
#include <vector>

namespace elephon
{
namespace LatticeStructure
{

template<class C>
class SymmetryReduction
{
public:
	SymmetryReduction();

	SymmetryReduction(
			LatticeStructure::Symmetry const & sym,
			std::vector<C> const & reducible,
			std::vector<C> & irreducible,
			std::vector<int> & redToIrred,
			std::vector<int> & symRedToIrred,
			std::vector< std::vector<int> > & irredToRed,
			std::vector< std::vector<int> > & symIrredToRed );

	void reduce(
			LatticeStructure::Symmetry const & sym,
			std::vector<C> const & reducible,
			std::vector<C> & irreducible,
			std::vector<int> & redToIrred,
			std::vector<int> & symRedToIrred,
			std::vector< std::vector<int> > & irredToRed,
			std::vector< std::vector<int> > & symIrredToRed ) const;
};

} /* namespace LatticeStructure */
} /* namespace elephon */

#include "LatticeStructure/SymmetryReduction.hpp"
#endif /* ELEPHON_LATTICESTRUCTURE_SYMMETRYREDUCTION_H_ */
