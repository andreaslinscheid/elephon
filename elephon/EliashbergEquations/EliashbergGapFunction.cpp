/*	This file EliashbergGapFunction.cpp is part of elephon.
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
 *  Created on: May 16, 2018
 *      Author: A. Linscheid
 */

#include "EliashbergEquations/EliashbergGapFunction.h"

namespace elephon {
namespace EliashbergEquations {

void
EliashbergGapFunction::initialize(int nMats, int nBands, EliashbergModule::EliashbergDataType value)
{
	MatsubaraBaseFermion::initialize(nMats, nBands,
			Auxillary::alignedvector::aligned_vector<EliashbergModule::EliashbergDataType>(nMats*nBands,value));
}

} /* namespace EliashbergEquations */
} /* namespace elephon */
