/*	This file ASSymmetry.cpp is part of elephon.
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
 *  Created on: Jan 2, 2018
 *      Author: A. Linscheid
 */

#include "AtomicSite/ASSymmetry.h"

namespace elephon
{
namespace AtomicSite
{


void
ASSymmetry::initialize(
		int lmax,
		std::vector<int> symmetryOperations,
		std::vector<double> fractionalTranslations,
		LatticeStructure::LatticeModule lattice,
		std::vector<LatticeStructure::Atom> const & atoms)
{

}

} /* namespace AtomicSite */
} /* namespace elephon */
