/*	This file Polynom.h is part of elephon.
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

#ifndef ELEPHON_AUXILLARY_CUBICPOLYNOMIAL_H_
#define ELEPHON_AUXILLARY_CUBICPOLYNOMIAL_H_

namespace elephon
{
namespace Auxillary
{

/**	A cubic polynomial for real numbers.
 *
 *  If TData is complex, this class will effectively be two polynomials for the real and imaginary part
 *  Independently. This is not a polynomial in the complex sense.
 *  TData*TMesh must be available and of type TData.
 *  TData and TMesh must be constructible by '0.0';
 *
 * 	The formula is according to Wikipedia [http://en.wikipedia.org/wiki/Spline_interpolation].
 * 	The template parameter T is float or double.
 * 	The right (larger) border is not part of the range of definition.
 * 	A polynomial obeys a comparison hirachie. We say one cubic polynomial is smaller than the other
 * 	by performing the same comparison for each the range of definition values.
 */
template < typename TMesh, typename TData >
class CubicPolynomial
{
public:
	/**
	 * Constructor that sets the polynomial.
	 *
	 * Internally calls CubicPolynomial::initialize. Please see there.
	 */
	CubicPolynomial(TMesh x1, TMesh x2, TData y1, TData y2, TData k1, TData k2);

	/**
	 * Initialize the Polynomial
	 *
	 * @param x1 The infinum of the range of definition.
	 * @param x2 The suppremum of the range of definition.
	 * @param y1 Data value at x1.
	 * @param y2 Data value at x2.
	 * @param k1 Derivative of data at x1.
	 * @param k2 Derivative of data at x2.
	 */
	void initialize(TMesh x1, TMesh x2, TData y1, TData y2, TData k1, TData k2);

	/**
	 * Evaluate the polynomial at x.
	 *
	 *  We allow evaluation outside the range of definition.
	 * @param x The position where to evaluate the polynomial.
	 * @param value The value of the polynomial at x.
	 */
	void evaluate(TMesh x, TData &value) const;

	/**
	 * Evaluate the derivative of the polynomial at x.
	 *
	 *  We allow evaluation outside the range of definition.
	 * @param x The position where to evaluate the polynomial.
	 * @param value The value of the derivative of the polynomial at x.
	 */
	void evaluate_derivative(TMesh x, TData &value) const;

	/**
	 * Evaluate the second derivative of the polynomial at x.
	 *
	 *  We allow evaluation outside the range of definition.
	 * @param x The position where to evaluate the polynomial.
	 * @param value The value of the second derivative of the polynomial at x.
	 */
	void evaluate_second_derivative(TMesh x, TData &value) const;

	/**
	 * @return The derivative of the data at range infinium.
	 */
	TData derivative_at_range_min() const;

	/**
	 * @return The derivative of the data at range suppremum.
	 */
	TData derivative_at_range_max() const;

	/**
	 * @return The left, inclusive boundary of the range of definition.
	 */
	TMesh get_range_min() const;

	/**
	 * @return The right, exclusive boundary of the range of definition.
	 */
	TMesh get_range_max() const;
private:

	TMesh rangeMin_ = TMesh(0.0); /// smallest value within the range of definition

	TMesh rangeMax_ = TMesh(0.0); /// smallest value larger than the range of definition

	TData y1_ = TData(0.0); ///data value at x1_

	TData y2_ = TData(0.0); ///data value at x2_

	TData a_ = TData(0.0); ///Derivative parameterization a = k1*( x2 - x1 ) - ( y2 - y1 )

	TData b_ = TData(0.0); ///Derivative parameterization b = -k2*( x2 - x1 ) + ( y2 - y1 )

	TMesh interval_length() const;
};

} /* namespace Auxillary */
} /* namespace elephon */

#include "Auxillary/CubicPolynomial.hpp"
#endif /* ELEPHON_AUXILLARY_CUBICPOLYNOMIAL_H_ */
