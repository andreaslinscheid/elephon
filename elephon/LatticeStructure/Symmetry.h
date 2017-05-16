/*	This file Symmetry.h is part of elephon.
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

#ifndef ELEPHON_LATTICESTRUCTURE_SYMMETRY_H_
#define ELEPHON_LATTICESTRUCTURE_SYMMETRY_H_

#include <vector>

namespace elephon
{
namespace LatticeStructure
{

class Symmetry
{
public:

	Symmetry();

	void initialize(
			double symmPrec,
			std::vector<int> symmetries,
			std::vector<double> fractionalTranslations);

	void apply(int isym, std::vector<double> & field) const;

	int get_index_inverse(int isym);

	double get_symmetry_prec() const;

	int get_num_symmetries() const;

	void symmetry_reduction( std::vector<int> const& indicesDropped);

private:

	double symmPrec_ = 0;

	int numSymmetries_ = 0;

	std::vector<int> inverseMap_;

	std::vector<int> symmetries_;

	std::vector<double> fractTrans_;
};

} /* namespace LatticeStructure */
} /* namespace elephon */

#endif /* ELEPHON_LATTICESTRUCTURE_SYMMETRY_H_ */
