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
#include "LatticeStructure/RegularSymmetricGrid.h"
#include "Auxillary/AlignedVector.h"
#include <vector>
#include <complex>
#include <memory>

namespace elephon
{
namespace PhononStructure
{

/**
 * Structure to store the data and provide methods related to the linear displacement potential due to a lattice perturbation.
 *
 * @todo This function is in some sense similar to ForceConstantMatrix.
 * 		 These two must be the first target of a refactor once the general logic is working ...
 */
class DisplacementPotential
{
public:

	void build(std::shared_ptr<const LatticeStructure::UnitCell> unitCell,
			std::shared_ptr<const LatticeStructure::UnitCell> superCell,
			std::shared_ptr<const std::vector<LatticeStructure::AtomDisplacement>> irredDispl,
			LatticeStructure::RegularBareGrid unitcellGrid,
			LatticeStructure::RegularBareGrid const & supercellGrid,
			std::vector<double> const & potentialUC,
			std::vector< std::vector<double> > potentialDispl,
			std::vector<int> coarseGrainGrid);

	void compute_dvscf_q(
			std::vector<double> const & qVect,
			Auxillary::alignedvector::DV const & modes,
			Auxillary::alignedvector::ZV const & dynamicalMatrices,
			std::vector<double> const & masses,
			std::vector<double> const & rVectors,
			Auxillary::alignedvector::CV & dvscf,
			std::vector<Auxillary::alignedvector::CV> & buffer,
			double freqCutoff = 0.0) const;

	int RVectorLayout(int iRx, int iRy, int iRz ) const;

	LatticeStructure::RegularBareGrid const & get_real_space_grid() const;

	/**
	 * Request a connection of arbitrary q points to the closest point in the coarse grained q-grid.
	 *
	 * If coarse graining is inactive, the code will simply return a map with each q point associated
	 * to its index in the array qVect.
	 *
	 * @param qVect						non-grid q vectors in the layout q0x, q0y, q0z, q1x ...
	 * @param qGridCorseGainToQIndex	a map that tells for each coarse grid vectors the indices
	 * 									in the list of q vectors qVect that are closest to it.
	 */
	void query_q(
			std::vector<double> const & qVect,
			std::map<LatticeStructure::RegularBareGrid::GridPoint, std::vector<int>> & qGridCorseGainToQIndex) const;

	int get_num_R() const;

	int get_num_modes() const;

	void write_dvscf(int atomIndex, int xi, std::string filename) const;

	void write_dvscf_q(std::vector<double> const & qVect,
			std::vector<int> modeIndices,
			Auxillary::alignedvector::DV const & modes,
			Auxillary::alignedvector::ZV const & dynamicalMatrices,
			std::vector<double> const & masses,
			std::string filename) const;
private:

	int numModes_ = 0;

	int nptsRealSpace_ = 0;

	LatticeStructure::RegularBareGrid unitCellGrid_;

	std::shared_ptr<const LatticeStructure::UnitCell> unitCell_;

	std::shared_ptr<LatticeStructure::RegularBareGrid> coarseGrainGrid_;

	///For each lattice vector, the linear displacement potential sorted
	///according 1) the mode index and 2) to the realspaceGrid_ in z-slowest-running order.
	std::vector<float> data_;

	std::vector<int> superCellDim_;

	std::vector<int> RVectorDim_;

	/// A buffer for the R vectors. These are implicitly defined by superCellDim_
	std::vector<double> RVectors_;

	std::vector<int> numRSpacePointsForLatticeVector_;

	void compute_rot_map(
			std::vector<double> const & shift,
			LatticeStructure::RegularBareGrid const & supercellGrid,
			LatticeStructure::Symmetry const & siteSymmetry,
			std::vector< std::vector<int> > & rotMap) const;

	int mem_layout(int ir, int mu, int iR ) const;

	/**
	 * Establish a connection between a point in the supercell and the corresponding one in the primitive cell.
	 *
	 * This means to create a table which maps a point in real space in the supercell
	 * into a point in the primitive cell plus a lattice vector
	 * With 'primitive' we mean the point in the reference frame of the infinite lattice where a given atom
	 * is always in the center. Thus even for the trivial supercell, this will not be an identity mapping.
	 *
	 * @param primitiveCellGrid		The grid which describes the real space of the primitive cell
	 * @param supercellGrid			The grid which describes the real space of the supercell
	 * @param atomsUC				A list of atoms in the primitive cell.
	 * @param rSuperCellToPrimitve	a vector that will be resized to hold 1) for each atom 2) a list of
	 * 								pairs where the first (int) references the consequtive real space index
	 * 								in the primitive cell and the second is a vector of 3 ints labeling
	 * 								the R vectors
	 */
	void build_supercell_to_primitive(
			LatticeStructure::RegularBareGrid const & primitiveCellGrid,
			LatticeStructure::RegularBareGrid const & supercellGrid,
			std::vector<LatticeStructure::Atom> const & atomsUC,
			std::vector< std::vector< std::pair<int,std::vector<int> > > > & rSuperCellToPrimitve) const;

	void clean_displacement_potential();

	void set_R_vectors();

	/**
	 * This method subtracts from the potentialDispl the periodically repeated one of the prestive unit cell.
	 */
	void constuct_potential_variation(
			std::vector<double> const & potentialUC,
			int nRSC,
			LatticeStructure::RegularBareGrid const & unitcellGrid,
			LatticeStructure::RegularBareGrid const & supercellGrid,
			std::vector< std::vector<double> > & potentialDispl);

	void symmetrize_periodic_dvscf_q(
			std::vector<double> const & qpoints,
			Auxillary::alignedvector::CV & dvscfData) const;

	void multiply_phase_qr(
			std::vector<double> const & qpoints,
			Auxillary::alignedvector::CV &dvscf_q) const;

};

} /* namespace PhononStructure */
} /* namespace elephon */

#endif /* ELEPHON_PHONONSTRUCTURE_DISPLACEMENTPOTENTIAL_H_ */
