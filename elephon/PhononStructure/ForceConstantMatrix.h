/*	This file ForceConstantMatrix.h is part of elephon.
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
 *  Created on: May 31, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_PHONONSTRUCTURE_FORCECONSTANTMATRIX_H_
#define ELEPHON_PHONONSTRUCTURE_FORCECONSTANTMATRIX_H_

#include "LatticeStructure/UnitCell.h"
#include "LatticeStructure/AtomDisplacement.h"
#include "LatticeStructure/Symmetry.h"
#include <vector>
#include <complex>
#include <map>

namespace elephon
{
namespace PhononStructure
{

class ForceConstantMatrix
{
public:

	void build(   LatticeStructure::UnitCell unitCell,
			LatticeStructure::UnitCell const & superCell,
			std::vector<LatticeStructure::AtomDisplacement> const & irredDispl,
			std::vector<std::vector<double>> forces);

	double operator() (int mu1, int mu2, int Rx, int Ry, int Rz) const;

	void fourier_transform_q(std::vector<double> const & qVect,
			std::vector<std::complex<double>> & data) const;

	int get_num_modes() const;

	int get_num_R() const;
private:

	int numModes_ = 0;

	std::vector<int> supercellDim_;

	//The R vector grid spans the cell from [1,1] including both borders in each direction
	std::vector<int> RVectorDim_;

	std::vector<int> RSupercellMultiplicity_;

	std::vector<double> data_;

	std::vector<double> tau_;

	void set_supercell_dim(
			LatticeStructure::UnitCell const & unitCell,
			LatticeStructure::UnitCell const & superCell);

	void determine_lattice_vector(LatticeStructure::AtomDisplacement const & displ,
			std::vector<LatticeStructure::Atom> const & superCell,
			std::vector<int> & R, int & mu);

	int RVectorLayout(int iRz, int iRy, int iRx ) const;

	int mem_layout(int iRz, int iRy, int iRx, int mu2, int mu1 ) const;

	int mem_layout( int ir, int mu2, int mu1 ) const;

	void shift_atoms_into_map( std::vector<LatticeStructure::Atom> const & atoms,
			LatticeStructure::Atom const & site,
			std::map<LatticeStructure::Atom,int> & shiftedLattice) const;

	void transform_map( std::map<LatticeStructure::Atom,int> & shiftedLattice,
			LatticeStructure::Symmetry::SymmetryOperation const & symOp,
			std::vector<int> & rotAtomsMap) const;

	void drift_clean_forces(
			std::vector<std::vector<double>> forces) const;
};

} /* namespace PhononStructure */
} /* namespace elephon */

#endif /* ELEPHON_PHONONSTRUCTURE_FORCECONSTANTMATRIX_H_ */
