/*	This file CubeSplineInterpolation.h is part of elephon.
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

#ifndef ELEPHON_ALGORITHMS_CUBESPLINEINTERPOLATION_H_
#define ELEPHON_ALGORITHMS_CUBESPLINEINTERPOLATION_H_

#include "Auxillary/CubicPolynomial.h"
#include <vector>
#include <limits>
#include <memory>

namespace elephon
{
namespace Algorithms
{
/**
 * A cubic spline that interpolates data.
 *
 * 	The formula is according to Wikipedia [http://en.wikipedia.org/wiki/Spline_interpolation]
 *
 * 	Very unfortunately, the BOOST spine requires equally spaced support points which is not
 * 	an option for us. We attempt to provide comparable usage characteristics.
 */
template<typename TMesh = double,
		typename TData = double>
class CubeSplineInterpolation
{
public:

	/**
	 * Empty constructor calls just CubeSplineInterpolation::clear().
	 */
	CubeSplineInterpolation();

	/**
	 * Initialize the spline.
	 *
	 * The formula is according to wikipedia [https://en.wikipedia.org/wiki/Spline_interpolation]
	 *
	 * @param dataBegin		Random access iterator pointing to the beginning of the range support function values f_i
	 * @param dataEnd		Random acces iterator pointing to the end of the range for support function values f_i
	 * @param meshBegin		Positions x_i where the function values are defined.
	 * @param meshEnd		Positions x_i where the function values are defined.
	 * @param intervalBegin		Beginning of the range of definition
	 * @param intervalEnd		End of the range of definition
	 * @param splineMatrixPtr	A shared pointer. If not null, but the vector type pointed to is empty,
	 * 							 the inverse spline matrix will place behind the pointer
	 * 							 for reuse. If not null and the data pointed has size not empty, it must be the
	 * 							 correct inverse spline matrix for that problem.
	 */
	template<class RandomAccessIteratorData, class RandomAccessIteratorMesh>
	void initialize(
			RandomAccessIteratorData dataBegin,
			RandomAccessIteratorData dataEnd,
			RandomAccessIteratorMesh meshBegin,
			RandomAccessIteratorMesh meshEnd,
			TMesh intervalBegin,
			TMesh intervalEnd,
			std::shared_ptr<Auxillary::alignedvector::aligned_vector<TMesh>> splineMatrixPtr = nullptr);

	/**
	 * Initialize the spline though the high level interface.
	 *
	 * @param mesh 		A vector type with the mesh data.
	 * @param data		A vector type with function values at the corresponding mesh point.
	 */
	template<class VTMesh, class VTData>
	void initialize(
			VTMesh const & mesh,
			VTData const & data);

	/**
	 * Removes the content of the object and sets it to its initial state.
	 */
	void clear();

	/**
	 * Evaluate the spline at point x.
	 *
	 * @param x		Location where to evaluate the spline.
	 * @return		the value of the spline at x
	 */
	TData operator() (TMesh x) const;

	/**
	 * Evaluate the derivative of the spline at point x.
	 *
	 * @param x		Location where to evaluate the derivative of the spline.
	 * @return		the value of the derivative of the spline at x.
	 */
	TData prime(TMesh x) const;

	/**
	 * @return The minimum of the defined range.
	 */
	TMesh min_range() const;

	/**
	 * @return The maximum of the defined range.
	 */
	TMesh max_range() const;
private:

	std::vector<Auxillary::CubicPolynomial<TMesh,TData>> polys_;

	TMesh rangeOfDefinitionBegin = 0.0;

	TMesh rangeOfDefinitionEnd = 0.0;

	int find_polynomial_index_in_range( TMesh x ) const;

	template<class RandomAccessIteratorMesh>
	void compute_inverse_spline_matrix(
			RandomAccessIteratorMesh meshBegin,
			RandomAccessIteratorMesh meshEnd,
			Auxillary::alignedvector::aligned_vector<TMesh> & inverseSplineMatrix ) const;
};

} /* namespace Algorithms */
} /* namespace elephon */

#include "Algorithms/CubeSplineInterpolation.hpp"
#endif /* ELEPHON_ALGORITHMS_CUBESPLINEINTERPOLATION_H_ */
