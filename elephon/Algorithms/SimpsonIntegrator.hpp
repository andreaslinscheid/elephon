/*	This file SimpsonIntegrator.hpp is part of elephon.
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
 *  Created on: May 17, 2018
 *      Author: A. Linscheid
 */

#include "Algorithms/SimpsonIntegrator.h"

namespace elephon {
namespace Algorithms {

template<class D>
typename SimpsonIntegrator<D>::T
SimpsonIntegrator<D>::integrate(typename D::const_iterator dataArrayBegin, typename D::const_iterator dataArrayEnd,
		typename D::const_iterator positionsBegin, typename D::const_iterator positionsEnd)
{
	assert(std::distance(positionsBegin, positionsEnd) == std::distance(dataArrayBegin, dataArrayEnd));
	if (std::distance(positionsBegin, positionsEnd) == 1)
		return T(0); // zero measure

	if (std::distance(positionsBegin, positionsEnd) == 2)
	{
		auto measure = -(*positionsBegin);
		measure += *(++positionsBegin);
		auto data = (*dataArrayBegin);
		data += *(++dataArrayBegin);
		return 0.5*measure*data; // Trapezoidal rule
	}

	T integral = T(0);
	auto p = positionsBegin;
	auto pm2 = p++;
	auto pm1 = p++;
	auto d = dataArrayBegin;
	auto dm2 = d++;
	auto dm1 = d++;
	for (;; ++ ++pm1, ++ ++pm2, ++ ++d, ++ ++dm1, ++ ++dm2)
	{
		const T h0 = *pm1-*pm2;
		const T h1 = *p-*pm1;
		const T hsum = h0 + h1;
		const T hprod = h0 * h1;
		const T h0divh1 = h0 / h1;
		integral += hsum/6.0*((*dm2)*(2-1.0/h0divh1) + (*dm1)*hsum*hsum/hprod + (*d)*(2-h0divh1));
		++p; // increase once
		if ( p == positionsEnd ){
			break;
		}

		++p; // increase the second time
		if ( p == positionsEnd ){
			++pm1;
			++pm2;
			++dm1;
			++dm2;
			// add the uneven end interval
			integral += 0.5*(*pm1-*pm2)*(*dm1-*dm2);
			break;
		}
	}
	return integral;
}

template<class D>
typename SimpsonIntegrator<D>::T
SimpsonIntegrator<D>::integrate(D const & dataArray, D const & positions)
{
	return this->integrate(dataArray.begin(), dataArray.end(), positions.begin(), positions.end());
}

} /* namespace Algorithms */
} /* namespace elephon */
