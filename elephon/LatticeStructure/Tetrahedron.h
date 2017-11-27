/*	This file Tetrahedron.h is part of elephon.
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
 *  Created on: Nov 2, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_LATTICESTRUCTURE_TETRAHEDRON_H_
#define ELEPHON_LATTICESTRUCTURE_TETRAHEDRON_H_

#include "LatticeStructure/ExtendedSymmetricGrid.h"
#include <vector>
#include <memory>

namespace elephon
{
namespace LatticeStructure
{

class Tetrahedron
{
public:
	Tetrahedron(
			std::vector<int> cornerIndicesData,
			std::vector<int> cornerIndicesExtended,
			std::vector<int> cornerIndicesExtendedReducible,
			std::shared_ptr<const ExtendedSymmetricGrid> extendedGrid);

	int get_multiplicity() const;

	void set_multiplicity(int m);

	std::vector<int> const & get_corner_indices() const;

	void compute_corner_vectors(
			std::vector<double> & v123 ) const;

	void compute_corner_points(
			std::vector<double> & p0123 ) const;

	void check_vectors_inside(
			std::vector<double> const & v,
			std::vector<bool> & inside,
			std::vector<double> & barycentricCoordinates) const;

	friend bool operator< (Tetrahedron const & t1, Tetrahedron const & t2);

private:

	int multiplicity_ = 1;

	std::vector<int> cornerIndicesData_;

	std::vector<int> cornerIndicesExtended_;

	std::vector<int> cornerIndicesExtendedReducible_;

	std::shared_ptr<const ExtendedSymmetricGrid> extendedGrid_;
};


} /* namespace LatticeStructure */
} /* namespace elephon */

#endif /* ELEPHON_LATTICESTRUCTURE_TETRAHEDRON_H_ */
