/*	This file RadialGrid.hpp is part of elephon.
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
 *  Created on: Jan 4, 2018
 *      Author: A. Linscheid
 */

#include "AtomicSite/RadialGrid.h"
#include "Auxillary/AlignedVector.h"
#include "Algorithms/CubeSplineInterpolation.h"
#include <cassert>
#include <iterator>
#include <memory>

namespace elephon
{
namespace AtomicSite
{

template<class ConstIterator, class Iterator>
void
RadialGrid::interpolate(
	std::vector<double> const & rValues,
	int nDataPerRValue,
	ConstIterator gridDataBegin,
	ConstIterator gridDataEnd,
	Iterator interpolDataBegin ) const
{
	const int nPInterpol = rValues.size();
	const int nPData = radialPoints_.size();
	assert(std::distance(gridDataBegin, gridDataEnd) == nPData);
	assert(numR_%nPData == 0);
	assert(nDataPerRValue == nPInterpol/nPData);

	// perform spline interpolation, we keep the spline matrix for efficiency
	auto splineMatrixPtr = std::make_shared<Auxillary::alignedvector::aligned_vector<double>>();
	for (int iBlock = 0 ; iBlock < nDataPerRValue; ++iBlock )
	{
		ConstIterator gridDataBeginBlock = gridDataBegin + iBlock*nPData;
		ConstIterator gridDataEndBlock = gridDataBegin + (iBlock+1)*nPData;
		Algorithms::CubeSplineInterpolation< double,
											 typename std::iterator_traits<ConstIterator>::value_type
											> csinterpol;
		csinterpol.initialize(
				gridDataBeginBlock, gridDataEndBlock,
				radialPoints_.begin(), radialPoints_.end(),
				0.0, radius_,
				splineMatrixPtr);

		for (auto r : rValues)
		{
			*(interpolDataBegin + nPInterpol*iBlock) = csinterpol(r);
		}
	}
}


} /* namespace AtomicSite */
} /* namespace elephon */
