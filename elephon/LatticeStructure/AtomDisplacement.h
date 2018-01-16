/*	This file AtomDisplacement.h is part of elephon.
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
 *  Created on: May 15, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_LATTICESTRUCTURE_ATOMDISPLACEMENT_H_
#define ELEPHON_LATTICESTRUCTURE_ATOMDISPLACEMENT_H_

#include "Symmetry.h"
#include "symmetry/SymmetryOperation.h"
#include <vector>
#include <string>

namespace elephon
{
namespace LatticeStructure
{

class AtomDisplacement
{
public:

	AtomDisplacement();

	AtomDisplacement(
			std::string kind,
			double magnitude,
			std::vector<double> position,
			std::vector<double> direction,
			double gridPrecision = 1e-6,
			bool symmetricDirection = true);

	void initialize(
			std::string kind,
			double magnitude,
			std::vector<double> position,
			std::vector<double> direction,
			double gridPrecision = 1e-6,
			bool symmetricDirection = true);

	void scale_position(double scaleX, double scaleY, double scaleZ);

	void transform(symmetry::SymmetryOperation const& sop);

	void transform_direction(symmetry::SymmetryOperation const& sop);

	double get_prec() const;

	std::vector<double> const & get_position() const;

	std::vector<double> const & get_direction() const;

	std::string get_kind() const;

	double get_magnitude() const;

	bool is_plus_minus_displ_equivalent() const;

	//Strict weak ordering is implied by (in that order)
	//	kind then position [x,y,z] then direction [|x|,|y|,|z|]
	//NOTE: direction is compared in magnitudes of the components.
	//		since the derivative is symmetric a displacement in direction (-1,0,0)
	//		is considered equal to (1,0,0).
	//Example: 	We compare displacements '1' and '2'.
	//			First we compare the kinds according to string.compare()
	//			Then position[x] is compared between to displacements '1' and '2'.
	//			If |A1 - A2| > equivalencePrc_ this implies '1' < '2'.
	//			If they are (almost) equal, we compare position y and so on.
	// This implies '1' == '2' if neither '1'<'2' nor '2' < '1'.
	friend bool operator< (AtomDisplacement const& d1, AtomDisplacement const& d2);


	friend bool operator== (AtomDisplacement const& d1, AtomDisplacement const& d2);

private:

	double equivalencePrc_ = 1e-6;

	bool treatDirectionSymmetric_ = true;

	std::string kind_;

	double magnitude_ = 0;

	std::vector<double> position_ = {0.0,0.0,0.0};

	//In cartesian coordinates
	std::vector<double> direction_ = {1.0,0.0,0.0};

	void normalize_direction();
};

} /* namespace LatticeStructure */
} /* namespace elephon */

#endif /* ELEPHON_LATTICESTRUCTURE_ATOMDISPLACEMENT_H_ */
