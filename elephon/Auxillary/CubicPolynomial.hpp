/*	This file CubicPolynomial.hpp is part of elephon.
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
 *  Created on: Jan 7, 2018
 *      Author: A. Linscheid
 */

#include "Auxillary/CubicPolynomial.h"
#include <cmath>
#include <cassert>

namespace elephon
{
namespace Auxillary
{

template<typename TMesh, typename TData>
CubicPolynomial<TMesh, TData>::CubicPolynomial(
		TMesh x1, TMesh x2,
		TData y1, TData y2,
		TData k1, TData k2)
{
	this->initialize(x1,x2,y1,y2,k1,k2);
}

template<typename TMesh, typename TData>
void
CubicPolynomial<TMesh, TData>::initialize(
		TMesh x1, TMesh x2,
		TData y1, TData y2,
		TData k1, TData k2)
{
	rangeMin_ = x1;
	rangeMax_ = x2;
	assert(rangeMin_ <= rangeMax_);
	a_ = k1*TData( x2 - x1 ) - ( y2 - y1 );
	b_ = -k2*TData( x2 - x1 ) + ( y2 - y1 );
	y1_ = y1;
	y2_ = y2;
}

template<typename TMesh, typename TData>
void
CubicPolynomial<TMesh, TData>::evaluate(
		TMesh x, TData &value) const
{
	TMesh t = ( x - rangeMin_ ) / this->interval_length() ;
	value = y1_*TData(1-t)+y2_*TData(t)+(a_*TData(1-t)+b_*TData(t))*TData(t*(1-t));
}

template<typename TMesh, typename TData>
void
CubicPolynomial<TMesh, TData>::evaluate_derivative(
		TMesh x, TData &value) const
{
	TMesh t = ( x - rangeMin_ ) / this->interval_length() ;
	value = (y2_-y1_)*TData(1.0/this->interval_length())+(a_*TData(1-t)+b_*TData(t))*TData(1-2*t)
			*TData(1.0/this->interval_length()) + (b_-a_)*TData(t*(1-t)/this->interval_length());
}

template<typename TMesh, typename TData>
void
CubicPolynomial<TMesh, TData>::evaluate_second_derivative(
		TMesh x, TData &value) const
{
	TMesh t = ( x - rangeMin_ ) / this->interval_length();
	value = (b_-2*a_+(a_-b_)*TData(3*t))*TData(2.0/std::pow(this->interval_length(),2));
}

template<typename TMesh, typename TData>
TData
CubicPolynomial<TMesh, TData>::derivative_at_range_min() const
{
	return (y2_-y1_+a_)*(1.0/this->interval_length());
}

template<typename TMesh, typename TData>
TData
CubicPolynomial<TMesh, TData>::derivative_at_range_max() const
{
	return (y2_-y1_-b_)/this->interval_length();
}

template<typename TMesh, typename TData>
TMesh
CubicPolynomial<TMesh, TData>::interval_length() const
{
	return rangeMax_-rangeMin_;
}

template<typename TMesh, typename TData>
TMesh
CubicPolynomial<TMesh, TData>::get_range_min() const
{
	return rangeMin_;
}

template<typename TMesh, typename TData>
TMesh
CubicPolynomial<TMesh, TData>::get_range_max() const
{
	return rangeMax_;
}

} /* namespace Auxillary */
} /* namespace elephon */
