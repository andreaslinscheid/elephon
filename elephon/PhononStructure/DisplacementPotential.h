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
#include "LatticeStructure/RegularSymmetricGrid.h"
#include "Auxillary/AlignedVector.h"
#include <vector>
#include <complex>
#include <memory>

namespace elephon
{
namespace LatticeStructure {class AtomDisplacementCollection; };
namespace LatticeStructure {class PrimitiveToSupercellConnection; };
namespace PhononStructure
{

class PotentialChangeIrredDisplacement;

/**
 * 	Structure to store the data and provide methods related to the linear displacement potential due to a lattice perturbation.
 *
 *	We denote the object this class handles as the row vector
 *	\f{eqnarray*}{
 *	{\bf{\delta v}_{{\rm {\scriptscriptstyle scf}}}^{\Sigma}}_{\kappa^{\prime}}(\bf{r}-\bf{R})
 *	& = &
 *	(\delta v_{{\rm {\scriptscriptstyle scf}}}^{\kappa x}(\bf{r}-\bf{R})\;\delta v_{{\rm {\scriptscriptstyle scf}}}^{\kappa y}(\bf{r}-\bf{R})\;
 *	\delta v_{{\rm {\scriptscriptstyle scf}}}^{\kappa z})(\bf{r}-\bf{R})
 *	\f}
 *	The central functionalities are first compute this data from irreducible input data in DisplacementPotential::initialize()
 *	and second compute the lattice Fourier DisplacementPotential::compute_dvscf_q().
 */
class DisplacementPotential
{
public:

	/**
	 * Compute \f$ {\bf{\delta v}_{{\rm {\scriptscriptstyle scf}}}^{\Sigma}}_{\kappa^{\prime}}(\bf{r}-\bf{R})\f$ from
	 * irreducible input data and store it in this object.
	 *
	 * @param[in] unitCell
	 * @param[in] superCell
	 * @param[in] displCollection
	 * @param[in] primToSCCon
	 * @param[in] unitcellGrid
	 * @param[in] supercellGrid
	 * @param[in] potentialChange
	 */
	void initialize(
			std::shared_ptr<const LatticeStructure::UnitCell> unitCell,
			std::shared_ptr<const LatticeStructure::UnitCell> superCell,
			std::shared_ptr<const LatticeStructure::AtomDisplacementCollection> displCollection,
			std::shared_ptr<const LatticeStructure::PrimitiveToSupercellConnection> primToSCCon,
			LatticeStructure::RegularBareGrid unitcellGrid,
			LatticeStructure::RegularBareGrid const & supercellGrid,
			std::vector<std::shared_ptr<const PotentialChangeIrredDisplacement>> const & potentialChange);

	/**
	 * Compute \f$ \delta v_{{\rm {\scriptscriptstyle scf}}}(\bf{r}) / \delta u_{\bf{q}\nu} \f$ from
	 * data stored in this object at a requested list of q vectors.
	 *
	 * @param[in] qVect
	 * @param[in] modes
	 * @param[in] dynamicalMatrices
	 * @param[in] masses
	 * @param[out] dvscf
	 * @param[in,out] buffer
	 * @param[in] freqCutoff
	 */
	void compute_dvscf_q(
			std::vector<double> const & qVect,
			Auxillary::Multi_array<double,2> const & modes,
			Auxillary::Multi_array<std::complex<double>,3> const & dynamicalMatrices,
			std::vector<double> const & masses,
			Auxillary::alignedvector::CV & dvscf,
			std::vector<Auxillary::alignedvector::CV> & buffer,
			double freqCutoff = 0.0) const;

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
			Auxillary::Multi_array<double,2> const & modes,
			Auxillary::Multi_array<std::complex<double>,3> const & dynamicalMatrices,
			std::vector<double> const & masses,
			std::string filename) const;

	Auxillary::Multi_array<float,3> const & get_data_regular_grid() const;

	Auxillary::Multi_array<std::complex<float>,4> const & get_data_radial_grid() const;
private:

	int numModes_ = 0;

	int nptsRealSpace_ = 0;

	LatticeStructure::RegularBareGrid primitiveCellGrid_;

	LatticeStructure::RegularBareGrid embeddingSuperCellGrid_;

	std::shared_ptr<const LatticeStructure::UnitCell> primitiveCell_;

	std::shared_ptr<const LatticeStructure::PrimitiveToSupercellConnection> primToSuperCell_;

	/// Regular grid part of the data. Coordinates are
	/// 0 : Lattice Vector index
	/// 1 : mode index, as layed out by Auxillary::memlayout::mode_layout()
	/// 2 : primitive Cell index
	Auxillary::Multi_array<float,3> dataRegular_;

	/// Radial grid part of the data. Coordinates are
	/// 0 : Lattice Vector index
	/// 1 : mode index;
	/// 2 : angular momentum + magnetic quantum number, as layed out by Auxillary::memlayout::angular_momentum_layout()
	/// 3 : radial data index;
	Auxillary::Multi_array<std::complex<float>,4> dataRadial_;


	/// Grid map of supercell indices to primtive cell indices, where the displaced atom is in the center. Coordinates are
	/// 0 : Atom index
	/// 1 : supercell grid index
	Auxillary::Multi_array<int, 2> perAtomGridIndexMap_;

	/// Grid map of supercell indices to lattice vectors to the cell from the primitive cell,
	/// where the displaced atom is in the center. Coordinates are
	/// 0 : Atom index
	/// 1 : supercell grid index
	/// 2 : R vector coordinates x y and z
	Auxillary::Multi_array<int, 3> perAtomGridIndexMapRVector_;

	std::vector<int> numRSpacePointsForLatticeVector_;

	/// For x,y and z the min and max of lattice vectors that occur in the embedding supercell.
	std::array<int, 3> embeddingRVectorBox_;

	/// A buffer for the embedding R vectors with the first coordinate is the index and the second of dim 3 has the x,y and z
	/// coordinates in units of the primitive cell. Note that in contrast to LatticeStructure::PrimitiveToSupercellConnection
	/// coordinates are such that the primitive cell is in the symmetric middle of the supercell.
	Auxillary::Multi_array<int,2> RVectors_;

	/// A mapping that tells for a tuple of Rx, Ry, and Rz coordintes within the ranges of embeddingRVectorBox_
	///  the index in the array RVectors_ or -1 if that index does not occur.
	Auxillary::Multi_array<int,3> RVectoToIndexMap_;

	int mem_layout(int ir, int mu, int iR ) const;

	void clean_displacement_potential();

	/**
	 * Given potential change data for irreducible displacements of an atom,
	 * generate the data for the full reducible set of displacements.
	 *
	 *
	 * @param siteSymmetry				The site symmetry operations of the atom
	 * @param potentialChange
	 * @param deltaVRegularSymExpanded
	 * @param deltaVRadialSymExpanded
	 */
	void symmetry_expand_displacement_data(
			LatticeStructure::Symmetry const & siteSymmetry,
			std::shared_ptr<const LatticeStructure::AtomDisplacementCollection> displColl,
			int atomIndex,
			int numAtomsSupercell,
			std::vector<std::shared_ptr<const PotentialChangeIrredDisplacement>> const & potentialChange,
			Auxillary::Multi_array<double,2> & deltaVRegularSymExpanded,
			Auxillary::Multi_array<std::complex<double>,4> & deltaVRadialSymExpanded) const;

	void symmetrize_periodic_dvscf_q(
			std::vector<double> const & qpoints,
			Auxillary::alignedvector::CV & dvscfData) const;

	void multiply_phase_qr(
			std::vector<double> const & qpoints,
			Auxillary::alignedvector::CV &dvscf_q) const;

	int R_vector_to_index(int Rx, int Ry, int Rz) const;

	void lattice_vector_map_back_embedding(std::array<int,3> & latticeVector) const;
};

} /* namespace PhononStructure */
} /* namespace elephon */

#endif /* ELEPHON_PHONONSTRUCTURE_DISPLACEMENTPOTENTIAL_H_ */
