/*	This file FrozenCore.h is part of elephon.
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
 *  Created on: Apr 27, 2018
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_ATOMICSITE_FROZENCORE_H_
#define ELEPHON_ATOMICSITE_FROZENCORE_H_

#include "AtomicSite/RadialGrid.h"
#include <vector>

namespace elephon
{
namespace AtomicSite
{

class FrozenCore
{
public:

	void initialize(double corePointCharge,
			std::vector<double> electronicFrozenCoreCharge,
			RadialGrid radialGrid);

private:

	double corePointCharge_ = 0.0;

	std::vector<double> electronicFrozenCoreCharge_;

	RadialGrid rgrid_;
};

} /* namespace AtomicSite */
} /* namespace elephon */

#endif /* ELEPHON_ATOMICSITE_FROZENCORE_H_ */
