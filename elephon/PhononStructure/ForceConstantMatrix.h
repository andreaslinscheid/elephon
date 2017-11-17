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
#include "Auxillary/AlignedVector.h"
#include <vector>
#include <complex>
#include <map>
#include <memory>

namespace elephon
{
namespace PhononStructure
{

class ForceConstantMatrix
{
public:

	void build( std::shared_ptr<const LatticeStructure::UnitCell> unitCell,
			std::shared_ptr<const LatticeStructure::UnitCell> superCell,
			std::shared_ptr<const std::vector<LatticeStructure::AtomDisplacement> > irredDispl,
			std::vector<std::vector<double>> forces);

	double operator() (int Rz, int Ry, int Rx,int mu2, int mu1) const;

	void fourier_transform_q(std::vector<double> const & qVect,
			Auxillary::alignedvector::ZV & data) const;

	void fourier_transform_derivative(std::vector<double> const & qVect,
			Auxillary::alignedvector::ZV & data) const;

	int get_num_modes() const;

	int get_num_R() const;

	void symmetrize(
			bool accusticSumRule,
			LatticeStructure::Symmetry const & symmetry );
private:

	int numModes_ = 0;

	std::vector<int> supercellDim_;

	std::vector<double> data_;

	std::shared_ptr<const LatticeStructure::UnitCell> uc_;

	///For all atoms 'a' in the unit cell, tau_ contains for all atoms 'b' in the unit cell
	// and the set of R vectors [(tau_a-tau_b)*numR+ir] a set with vectors of atom positions
	std::vector< std::vector<double> > tau_;

	void determine_lattice_vector(LatticeStructure::AtomDisplacement const & displ,
			std::vector<LatticeStructure::Atom> const & superCell,
			std::vector<int> & R, int & mu);

	int RVectorLayout(int iRd, int iRy, int iRz ) const;

	int mem_layout(int iRx, int iRy, int iRz, int mu2, int mu1 ) const;

	int mem_layout( int ir, int mu2, int mu1 ) const;

	void transform_map(
			std::vector<double> const & shift,
			std::vector<LatticeStructure::Atom> atoms,
			LatticeStructure::Symmetry const & siteSymmetry,
			std::vector< std::vector<int> > & rotAtomsMap) const;

	void drift_clean_forces(
			std::vector<std::vector<double>> forces) const;
};

} /* namespace PhononStructure */
} /* namespace elephon */

#endif /* ELEPHON_PHONONSTRUCTURE_FORCECONSTANTMATRIX_H_ */
