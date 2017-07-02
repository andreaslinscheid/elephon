/*	This file DisplacementPotential.h is part of elephon.
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
 *  Created on: Jun 30, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_PHONONSTRUCTURE_DISPLACEMENTPOTENTIAL_H_
#define ELEPHON_PHONONSTRUCTURE_DISPLACEMENTPOTENTIAL_H_

#include "LatticeStructure/UnitCell.h"
#include "LatticeStructure/AtomDisplacement.h"
#include "LatticeStructure/RegularGrid.h"

namespace elephon
{
namespace PhononStructure
{

class DisplacementPotential
{
public:

	void build(  LatticeStructure::UnitCell const & unitCell,
			LatticeStructure::UnitCell const & superCell,
			std::vector<LatticeStructure::AtomDisplacement> const & irredDispl,
			LatticeStructure::RegularGrid const & unitcellGrid,
			LatticeStructure::RegularGrid const & supercellGrid,
			std::vector<double> const & potentialUC,
			std::vector< std::vector<double> > const & potentialDispl );

	void compute_dvscf(
			std::vector<double> const & qVect,
			std::vector<double> const & dynamicalMatrices,
			std::vector<double> const & masses,
			std::vector<float> & dvscf) const;

	int RVectorLayout(int iRz, int iRy, int iRx ) const;

	int get_num_R() const;

	int get_num_modes() const;
private:

	int numModes_ = 0;

	int nptsRealSpace_;

	///For each lattice vector, the linear displacement potential sorted
	///according 1) the mode index and 2) to the realspaceGrid_ in z-slowest-running order.
	std::vector<float> data_;

	std::vector<int> superCellDim_;

	void set_R_vector_map(
			LatticeStructure::UnitCell const & unitCell,
			LatticeStructure::UnitCell const & superCell);

	void compute_rot_map(
			std::vector<double> const & shift,
			LatticeStructure::RegularGrid const & supercellGrid,
			LatticeStructure::Symmetry const & siteSymmetry,
			std::vector< std::vector<int> > & rotMap) const;

	int mem_layout(int ir, int mu, int iR ) const;
};

} /* namespace PhononStructure */
} /* namespace elephon */

#endif /* ELEPHON_PHONONSTRUCTURE_DISPLACEMENTPOTENTIAL_H_ */
