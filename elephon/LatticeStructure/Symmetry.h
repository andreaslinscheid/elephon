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

#include "LatticeModule.h"
#include <vector>
#include <assert.h>

namespace elephon
{
namespace LatticeStructure
{

class Symmetry
{
public:

	typedef struct Sop
	{
		int ptgroup[9];
		double fracTrans[3];
		void apply( std::vector<double> & v, bool latticePeriodic = true) const;
	} SymmetryOperation;

	Symmetry();

	void initialize(
			double symmPrec,
			std::vector<int> symmetries,
			std::vector<double> fractionalTranslations,
			LatticeStructure::LatticeModule lattice,
			bool hasTimeReversal);

	void set_reciprocal_space_sym();

	void apply(int isym, std::vector<double> & field, bool latticePeriodic = true) const;

	int get_index_inverse(int isym) const;

	int get_group_product(int isym1, int isym2) const;

	double get_symmetry_prec() const;

	bool is_symmorphic(int isym) const;

	bool has_inversion() const;

	bool is_inversion(int isym) const;

	int get_num_symmetries() const;

	void symmetry_reduction( std::vector<int> indicesDropped);

	SymmetryOperation get_sym_op( int isym ) const;

	int get_identity_index() const;

	bool is_reci() const;

	std::vector<double> get_fractional_translation(int isym) const;
private:

	int idIndex_ = 0;

	bool isReciprocalSpace_ = false;

	double symmPrec_ = 0;

	int numSymmetries_ = 0;

	int numRotations_ = 0;

	bool hasTimeReversal_ = true;

	bool hasInversion_ = true;

	std::vector<bool> isSymmorphic_;

	std::vector<int> inverseMap_;

	std::vector<int> multiplicationTable_;

	std::vector<int> symmetries_;

	std::vector<double> fractTrans_;

	std::vector<double> fractTransStore_;

	LatticeStructure::LatticeModule lattice_;
};

} /* namespace LatticeStructure */
} /* namespace elephon */

#endif /* ELEPHON_LATTICESTRUCTURE_SYMMETRY_H_ */
