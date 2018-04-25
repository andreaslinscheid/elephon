/*	This file ForceConstantMatrix.cpp is part of elephon.
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

#include "ForceConstantMatrix.h"
#include "Algorithms/LinearAlgebraInterface.h"
#include "LatticeStructure/Atom.h"
#include "LatticeStructure/SymmetryReduction.h"
#include "symmetry/atom_transform_map.h"
#include "LatticeStructure/UnitCell.h"
#include "LatticeStructure/AtomSymmetryConnection.h"
#include "symmetry/atom_transform_map.h"
#include "LatticeStructure/UnitCell.h"
#include "LatticeStructure/AtomDisplacementCollection.h"
#include "LatticeStructure/Symmetry.h"
#include "PhononStructure/Forces.h"
#include "LatticeStructure/PrimitiveToSupercellConnection.h"
#include "LatticeStructure/UnitCell.h"
#include <assert.h>
#include <cmath>
#include <set>
#include <map>
#include <iostream>
#include <memory>

namespace elephon
{
namespace PhononStructure
{

void
ForceConstantMatrix::initialize(std::shared_ptr<const LatticeStructure::UnitCell> primitiveCell,
		std::shared_ptr<const LatticeStructure::UnitCell> superCell,
		std::shared_ptr<const LatticeStructure::AtomDisplacementCollection> irredDispl,
		std::shared_ptr<const LatticeStructure::PrimitiveToSupercellConnection> primToSupercell,
		std::shared_ptr<const PhononStructure::Forces> forces)
{
	assert( forces->get_num_total_irred_displacements() == irredDispl->get_tota_num_irred_displacements() );
	primitiveCell_ = primitiveCell;
	this->set_tau_vectors_and_multiplicity(
			primitiveCell,
			superCell,
			primToSupercell);
	const int nR = this->get_num_R();
	const int nASC = superCell->get_atoms_list().size();
	const int nAPC = primitiveCell->get_atoms_list().size();
	numModes_ = nAPC*3;

	Algorithms::LinearAlgebraInterface linAlg;
	Auxillary::Multi_array<double,2> matrixForceConstants(boost::extents[3][3]);

	// We compute the matrix of force constants using the (L)ocal (M)atrix (L)ayout
	// before we reshape in the end to a format that is faster to be Fourier transformed
	// in the space of modes.
	Auxillary::Multi_array<double,5> dataLML(boost::extents[nR][nAPC][nAPC][3][3]);
	for ( int atomIndex : primitiveCell->get_atom_symmetry()->get_list_irreducible_atoms() )
	{
		LatticeStructure::Symmetry const & siteSym = primitiveCell->get_site_symmetry(atomIndex);

		// copy the irreducible displacement data for this atom
		int nRedDispl = irredDispl->get_num_red_displacements_for_atom(atomIndex);

		// symmetry - expand the irreducible displacements and the forces
		Auxillary::Multi_array<double,2> uInv;
		irredDispl->generate_pseudo_inverse_reducible_displacements_for_atom(atomIndex,uInv);
		Auxillary::Multi_array<double,3> forceSymExpanded;
		const int aSC = primToSupercell->primitive_to_supercell_atom_index(atomIndex);
		forces->site_symmetry_expand_data(siteSym, atomIndex, aSC, superCell->get_atoms_list(), forceSymExpanded);

		// multiplying the Moore-Penrose pseudo inverse, of the displacements,
		// we obtain the least square fit of transpose of the force constant matrix
		Auxillary::Multi_array<double,2> forcesThisAtom(boost::extents[nRedDispl][3]);
		for ( int iaS = 0 ; iaS < nASC ; ++iaS )
		{
			// We have to be a bit careful with the data layout. The forces as loaded have
			// for each displacement, for each atom in the supercell the forces in x, y and z
			// Here we take for a given atom in the supercell for each displacement the force in x, y and z
			// and multiply the inverse displacement matrix which gives our 3x3 matrix of force constants at this atom.
			for (int iRedDispl = 0 ; iRedDispl < nRedDispl; ++iRedDispl)
				for (int xi = 0 ; xi < 3; ++xi)
					forcesThisAtom[iRedDispl][xi] = forceSymExpanded[iRedDispl][iaS][xi];

			linAlg.matrix_matrix_prod( uInv, forcesThisAtom, matrixForceConstants, 3, 3 );

			const int iAPCEquiv = primToSupercell->get_equiv_atom_primitive(iaS);
			const int iR = primToSupercell->get_equiv_atom_primitive_lattice_vector_index(iaS);
			//transpose back
			for ( int xi = 0 ; xi < 3; ++xi)
				for ( int xj = 0 ; xj < 3; ++xj)
					dataLML[iR][atomIndex][iAPCEquiv][xi][xj] = -matrixForceConstants[xj][xi];
		}

		// transform the data to the star of symmetry equivalent atoms in the primitive cell.
		for ( auto atomStar : primitiveCell->get_atom_symmetry()->get_star_atom_indices(atomIndex) )
		{
			const int atomIndex_rot = atomStar.first;
			const int rSymmetryIndexIrredToRed = atomStar.second;
			if ( rSymmetryIndexIrredToRed != primitiveCell->get_symmetry().get_identity_index())
			{
				for ( int iaS = 0 ; iaS < nASC ; ++iaS )
				{
					const int iR = primToSupercell->get_equiv_atom_primitive_lattice_vector_index(iaS);
					const int iAPCEquiv = primToSupercell->get_equiv_atom_primitive(iaS);
					const int rotatedSCAtomIndex = superCell->get_atom_symmetry()->atom_rot_map(rSymmetryIndexIrredToRed, iaS);
					const int iR_rot = primToSupercell->get_equiv_atom_primitive_lattice_vector_index(rotatedSCAtomIndex);
					const int iAPCEquiv_rot = primToSupercell->get_equiv_atom_primitive(rotatedSCAtomIndex);

					for ( int xi = 0 ; xi < 3; ++xi)
						for ( int xj = 0 ; xj < 3; ++xj)
							matrixForceConstants[xi][xj] = dataLML[iR][atomIndex][iAPCEquiv][xi][xj];

					primitiveCell->get_symmetry().rotate_matrix_cartesian(
							rSymmetryIndexIrredToRed,
							matrixForceConstants.data(), matrixForceConstants.data()+matrixForceConstants.size());

					for ( int xi = 0 ; xi < 3; ++xi)
						for ( int xj = 0 ; xj < 3; ++xj)
							dataLML[iR_rot][atomIndex_rot][iAPCEquiv_rot][xi][xj] = matrixForceConstants[xi][xj];
				}
			}
		}
	}
	data_.resize(boost::extents[nR][numModes_][numModes_]);

	for (int iR = 0 ; iR < nR ; ++iR)
		for (int ia = 0; ia < nAPC; ++ia)
			for (int ib = 0; ib < nAPC; ++ib)
				for ( int xi = 0 ; xi < 3; ++xi)
					for ( int xj = 0 ; xj < 3; ++xj)
					{
						int mu1 = Auxillary::memlayout::mode_layout(ia,xi);
						int mu2 = Auxillary::memlayout::mode_layout(ib,xj);
						data_[iR][mu1][mu2] = dataLML[iR][ia][ib][xi][xj];
					}
}

double
ForceConstantMatrix::operator() (int Rx, int Ry, int Rz, int mu1, int mu2) const
{
	const int iR = indexLatticVectorMap_[Rx][Ry][Rz];
	assert((iR>=0)&&(iR<data_.shape()[0]));
	assert(((mu1>=0)&&(mu1<numModes_)) && ((mu2>=0)&&(mu2<numModes_)));
	return data_[iR][mu1][mu2];
}

int
ForceConstantMatrix::get_num_modes() const
{
	return numModes_;
}

int
ForceConstantMatrix::get_num_R() const
{
	return tau_.shape()[2];
}

void
ForceConstantMatrix::fourier_transform_q(
		std::vector<double> const & qVect,
		Auxillary::Multi_array<std::complex<double>,3> & dataFT,
		bool symmetrize) const
{
	assert( qVect.size()%3 == 0 );
	int nq = int(qVect.size())/3;
	int nM = this->get_num_modes();

	dataFT.resize(boost::extents[nq][nM][nM]);
	std::fill_n(dataFT.data(), dataFT.size(), std::complex<double>(0));

	int nAUC = nM/3;
	int nR = this->get_num_R();

	Auxillary::Multi_array<std::complex<double>,3> phases(boost::extents[nAUC][nAUC][nR]);
	for ( int iq = 0 ; iq < nq ; ++iq)
	{
		std::fill_n(phases.data(), phases.size(), std::complex<double>(0));
		for ( int ia = 0 ; ia < nAUC; ++ia )
			for ( int ib = 0 ; ib < nAUC; ++ib )
				for ( int iR = 0 ; iR < nR; ++iR )
					for ( int im = 0 ; im < multiplicity_[ia][ib][iR]; ++im )
						phases[ia][ib][iR] +=	std::exp( std::complex<double>(0,
										-2*M_PI*(qVect[iq*3+0]*tau_[ia][ib][iR][im][0]
											    +qVect[iq*3+1]*tau_[ia][ib][iR][im][1]
												+qVect[iq*3+2]*tau_[ia][ib][iR][im][2])) )
										/ double(multiplicity_[ia][ib][iR]);

		//only compute upper half + diagonal
		for ( int mu1 = 0 ; mu1 < nM; ++mu1 )
			for ( int mu2 = mu1 ; mu2 < nM; ++mu2 )
			{
				int ia = Auxillary::memlayout::atomIndex_of_mode(mu1);
				int ib = Auxillary::memlayout::atomIndex_of_mode(mu2);
				for ( int iR = 0 ; iR < nR; ++iR )
					dataFT[iq][mu1][mu2] += phases[ia][ib][iR]*data_[iR][mu1][mu2];
			}

		//reconstruct lower half as the conjugate
		for ( int mu1 = 0 ; mu1 < nM; ++mu1 )
			for ( int mu2 = mu1+1 ; mu2 < nM; ++mu2 )
				dataFT[iq][mu2][mu1] = std::conj(dataFT[iq][mu1][mu2]);
	}

	//symmetrize if requested
	if ( symmetrize )
	{
		this->symmetrize_q(qVect, dataFT);
	}
}

void
ForceConstantMatrix::symmetrize_q(
		std::vector<double> const & qpoints,
		Auxillary::Multi_array<std::complex<double>,3> & dataFT,
		std::shared_ptr<const LatticeStructure::UnitCell> unitCell) const
{
	if (! unitCell)
		unitCell = primitiveCell_;

	const int nq = qpoints.size()/3;
	const int nM = this->get_num_modes();
	const int nA = nM/3;
	assert(dataFT.size() == nq*nM*nM);

	auto fullSymmetryGroup = unitCell->get_symmetry();

	auto fullSymmetryGroupReci = fullSymmetryGroup;
	fullSymmetryGroupReci.set_reciprocal_space_sym(true);
	Auxillary::alignedvector::ZV rotBuffer(3*3);
	for ( int iq = 0 ; iq < nq ; ++iq)
	{
		std::vector<double> q(&qpoints[iq*3], &qpoints[iq*3]+3);
		auto dropedSymmetries = fullSymmetryGroupReci.get_list_incompatible_symops_in_group(q);
		const int nsymQ = fullSymmetryGroupReci.get_num_symmetries_no_T_rev() - static_cast<int>(dropedSymmetries.size());
		assert((nsymQ>0) && (nsymQ <= fullSymmetryGroupReci.get_num_symmetries_no_T_rev()));

		// for the pure identity symmetry there is nothing to be done.
		if (nsymQ == 1)
			continue;

		std::vector<bool> visited(nA*nA, false);

		std::set<int> dropSym(dropedSymmetries.begin(), dropedSymmetries.end());

		std::vector<std::complex<double>> rotPhases(nsymQ);

		for (int ia = 0 ; ia < nA ; ++ia)
			for (int ib = 0 ; ib < nA ; ++ib)
			{
				if ( visited[ia*nA+ib] )
					continue;

				std::fill(rotBuffer.begin(), rotBuffer.end(), std::complex<double>(0.0));

				// Note: 	we use an inverse logic here with 'dropped', instead of 'kept', symmetries is
				//			because the initial atom mapping is done for indices of full symmetry group.
				//			This way, we loop until we find a symmetry that is still present.
				int isymq = 0;
				for (int isym = 0 ; isym < fullSymmetryGroup.get_num_symmetries(); ++isym)
				{
					if ( dropSym.find(isym) != dropSym.end())
						continue;

					assert(isymq < nsymQ);
					auto symOp = fullSymmetryGroup.get_sym_op(isym);

					int mapped_ia = unitCell->get_atom_symmetry()->atom_rot_map(isym, ia);
					int mapped_ib = unitCell->get_atom_symmetry()->atom_rot_map(isym, ib);

					double arg = 0;
					for (int xi = 0 ; xi < 3 ; ++xi)
						arg +=  2.0*M_PI*q[xi]*(unitCell->get_atoms_list()[mapped_ia].get_position()[xi]
									- unitCell->get_atoms_list()[mapped_ib].get_position()[xi]);
					rotPhases[isymq] = std::complex<double>(std::cos(arg), std::sin(arg));

					for ( int i = 0; i < 3; ++i)
						for ( int j = 0; j < 3; ++j)
							for ( int k = 0; k < 3; ++k)
								for ( int l = 0; l < 3; ++l)
								{
									int mu1 = Auxillary::memlayout::mode_layout(mapped_ia,k);
									int mu2 = Auxillary::memlayout::mode_layout(mapped_ib,l);
									rotBuffer[i*3+j] += rotPhases[isymq]*
												symOp.get_carth_rot_matrix(i,k)*dataFT[iq][mu1][mu2]*symOp.get_carth_rot_matrix(j,l);
								}
					++isymq;
				}
				assert(isymq == nsymQ);

				// fill the symmetrized matrix in all symmetry equivalent places
				isymq = 0;
				for (int isym = 0 ; isym < fullSymmetryGroup.get_num_symmetries(); ++isym)
				{
					if ( dropSym.find(isym) != dropSym.end())
						continue;
					assert(isymq < nsymQ);
					auto symOp = fullSymmetryGroup.get_sym_op(isym);

					int mapped_ia = unitCell->get_atom_symmetry()->atom_rot_map(isym, ia);
					int mapped_ib = unitCell->get_atom_symmetry()->atom_rot_map(isym, ib);

					for ( int i = 0; i < 3; ++i)
						for ( int j = 0; j < 3; ++j)
						{
							int mu1 = Auxillary::memlayout::mode_layout(mapped_ia,i);
							int mu2 = Auxillary::memlayout::mode_layout(mapped_ib,j);
							dataFT[iq][mu1][mu2] = 0;
							for ( int k = 0; k < 3; ++k)
								for ( int l = 0; l < 3; ++l)
								{
									// note: symOp.ptgCart is transposed which is the inverse matrix
									dataFT[iq][mu1][mu2] += std::conj(rotPhases[isymq]) / double(nsymQ)*
												symOp.get_carth_rot_matrix(k,i)*rotBuffer[k*3+l]*symOp.get_carth_rot_matrix(l,j);
								}
						}

					visited[mapped_ia*nA+mapped_ib] = true;
					++isymq;
				}
				assert(isymq == nsymQ);
			}
	}
}

void
ForceConstantMatrix::fourier_transform_derivative(std::vector<double> const & qVect,
		Auxillary::alignedvector::ZV & data) const
{
	assert( qVect.size()%3 == 0 );
	const int nR = this->get_num_R();
	const int nq = int(qVect.size())/3;

	data.assign( nq*numModes_*numModes_*3 , 0.0 );
	for ( int iq = 0 ; iq < nq ; ++iq)
	{
		for ( int i = 0 ; i < 3 ; ++i )
			for ( int mu1 = 0 ; mu1 < numModes_; ++mu1 )
				for ( int mu2 = 0 ; mu2 < numModes_; ++mu2 )
				{
					int ia = Auxillary::memlayout::atomIndex_of_mode(mu1);
					int ib = Auxillary::memlayout::atomIndex_of_mode(mu2);
					for ( int ir = 0 ; ir < nR; ++ir )
					{
						for ( int im = 0 ; im < multiplicity_[ia][ib][ir]; ++im )
						{
							std::complex<double> phase = std::exp( std::complex<double>(0,
													-2*M_PI*(qVect[iq*3+0]*tau_[ia][ib][ir][im][0]
															+qVect[iq*3+1]*tau_[ia][ib][ir][im][1]
															+qVect[iq*3+2]*tau_[ia][ib][ir][im][2])) )
													/ double(multiplicity_[ia][ib][ir]);

							data[((iq*3+i)*numModes_ + mu1)*numModes_+mu2] +=
										std::complex<double>(0,-tauCart_[ia][ib][ir][im][i])*
										phase* data_[ir][mu1][mu2];
						}
					}
				}
	}
}

void
ForceConstantMatrix::drift_clean_forces(
		std::vector<std::vector<double>> forces) const
{
	//Remove residual forces by subtracting the ones in the unperturbed cell
	//and removing the sum of all forces, which must be zero in theory
	int nIrrdDispl = int(forces.size())-1;
	for ( int irredDispl = 0 ;  irredDispl < nIrrdDispl; ++irredDispl )
	{
		double sumForces[3] = {0,0,0};
		for ( int i = 0 ; i < forces[irredDispl].size(); ++i)
		{
			forces[irredDispl][i] -= forces[nIrrdDispl][i];
			sumForces[i%3] += forces[irredDispl][i];
		}
		for ( int i = 0 ; i < forces[irredDispl].size(); ++i)
			forces[irredDispl][i] -= sumForces[i%3];
	}
}

void
ForceConstantMatrix::set_tau_vectors_and_multiplicity(
		std::shared_ptr<const LatticeStructure::UnitCell> primitiveCell,
		std::shared_ptr<const LatticeStructure::UnitCell> superCell,
		std::shared_ptr<const LatticeStructure::PrimitiveToSupercellConnection> primToSupercell )
{
	Auxillary::Multi_array<int,2> RVectors;
	primToSupercell->get_supercell_vectors(RVectors, indexLatticVectorMap_);
	const int nR = RVectors.shape()[0];
	const int nAPC = primitiveCell->get_atoms_list().size();
	const int nASC = superCell->get_atoms_list().size();
	tau_.resize( boost::extents[nAPC][nAPC][nR] );
	multiplicity_.resize(boost::extents[nAPC][nAPC][nR]);
	for ( int ia = 0; ia < nAPC; ++ia )
	{
		// shift all the atoms in the supercell such that the
		// atom, equivalent to ia, is in the center.
		auto atomsSC = superCell->get_atoms_list();
		const int iaSCEquivalent = primToSupercell->primitive_to_supercell_atom_index(ia);
		auto shift = superCell->get_atoms_list()[iaSCEquivalent].get_position();
		for ( auto &a : atomsSC)
		{
			auto pos = a.get_position();
			for ( int i = 0 ; i < 3 ; ++i )
				pos[i] -= shift[i];
			a.set_position(pos);
		}

		//iff an atoms is close to the boarder we have to place a copy
		//on the respective other one to preserve the local symmetry [see PRL]
		for ( int iASC = 0 ; iASC < nASC; ++iASC)
		{
			const int ib = primToSupercell->get_equiv_atom_primitive(iASC);
			const int iR = primToSupercell->get_equiv_atom_primitive_lattice_vector_index(iASC);
			auto pos = atomsSC[iASC].get_position();
			std::vector<std::vector<double>> vects(1,pos);
			for ( int xi = 0 ; xi < 3; ++xi)
				if ( std::abs(std::abs(pos[xi])-0.5) < 1e-4 )
				{
					auto copyVects = vects;
					for ( auto taus : copyVects )
					{
						taus[xi] *= -1;
						vects.push_back(taus);
					}
				}

			//copy and scale to the unit of the primitive cell
			multiplicity_[ia][ib][iR] = vects.size();
			tau_[ia][ib][iR].resize(boost::extents[vects.size()][3]);
			for ( int iM = 0 ; iM < vects.size(); ++iM )
			{
				primToSupercell->supercell_to_primitive_coordinates(vects[iM]);
				for ( int xi = 0 ; xi < 3; ++xi)
					tau_[ia][ib][iR][iM][xi] = vects[iM][xi];
			}
		}
	}

	tauCart_.resize(boost::extents[nAPC][nAPC][nR]);
	for (int ia = 0 ; ia < nAPC; ++ia)
		for (int ib = 0 ; ib < nAPC; ++ib)
			for ( int iR = 0 ; iR < nR; ++iR )
			{
				tauCart_[ia][ib][iR].resize(boost::extents[tau_[ia][ib][iR].shape()[0]][3]);
				std::copy_n(tau_[ia][ib][iR].data(), tau_[ia][ib][iR].size(), tauCart_[ia][ib][iR].data());
				primitiveCell_->get_lattice().direct_to_cartesian_angstroem( tauCart_[ia][ib][iR].data(), tauCart_[ia][ib][iR].size());
			}
}

} /* namespace PhononStructure */
} /* namespace elephon */
