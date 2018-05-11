/*	This file DisplacementPotential.cpp is part of elephon.
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

#include "PhononStructure/DisplacementPotential.h"
#include "PhononStructure/PotentialChangeIrredDisplacement.h"
#include "LatticeStructure/SymmetryReduction.h"
#include "LatticeStructure/Atom.h"
#include "LatticeStructure/AtomDisplacementCollection.h"
#include "LatticeStructure/PrimitiveToSupercellConnection.h"
#include "LatticeStructure/AtomSymmetryConnection.h"
#include "Algorithms/LinearAlgebraInterface.h"
#include "IOMethods/WriteVASPRealSpaceData.h"
#include "Algorithms/GridRotationMap.h"
#include "Auxillary/UnitConversion.h"
#include "AtomicSite/SphericalHarmonicExpansion.h"
#include <iostream>

#include <iomanip>

namespace elephon
{
namespace PhononStructure
{

void
DisplacementPotential::initialize(
		std::shared_ptr<const LatticeStructure::UnitCell> primCell,
		std::shared_ptr<const LatticeStructure::UnitCell> superCell,
		std::shared_ptr<const LatticeStructure::AtomDisplacementCollection> displCollection,
		std::shared_ptr<const LatticeStructure::PrimitiveToSupercellConnection> primToSCCon,
		LatticeStructure::RegularBareGrid primcellGrid,
		LatticeStructure::RegularBareGrid const & supercellGrid,
		std::vector<std::shared_ptr<const PotentialChangeIrredDisplacement>> const & potentialChange )
{
	assert(potentialChange.size()>0);

	primToSCCon->build_embedding_supercell(
			primcellGrid,
			supercellGrid,
			embeddingRVectorBox_,
			RVectors_,
			RVectoToIndexMap_,
			embeddingSuperCellGrid_);

	// precompute the supercell to primitive cell folding maps
	primToSCCon->build_supercell_to_primitive(
			primcellGrid,
			supercellGrid,
			primCell->get_atoms_list(),
			embeddingRVectorBox_,
			perAtomGridIndexMap_,
			perAtomGridIndexMapRVector_);

	numModes_ = primCell->get_atoms_list().size()*3;
	const int nrSC = supercellGrid.get_num_points();
	const int nrPC = primcellGrid.get_num_points();
	const int nASC = superCell->get_atoms_list().size();
	const int nR = this->get_num_R();
	const int nRadMax = potentialChange[0]->get_max_num_radial_elements();
	const int nAngChnlMax = potentialChange[0]->get_max_num_angular_moment_channels();
	const int LMax = potentialChange[0]->get_max_angular_moment();

	nptsRealSpace_ = nrPC;
	assert( primcellGrid.get_num_points() == nrPC );
	assert( potentialChange.size() == displCollection->get_tota_num_irred_displacements() );

	Algorithms::LinearAlgebraInterface linAlg;
	Auxillary::Multi_array<double,2> pseudoInverseDisplacements;
	Auxillary::Multi_array<std::complex<double>,2> pseudoInverseDisplacementsCmplx;
	Auxillary::Multi_array<double,2> deltaVRegularSymExpanded;
	Auxillary::Multi_array<std::complex<double>,4> deltaVRadialSymExpanded;
	Auxillary::alignedvector::DV linDVscfRegular, linDVscfRegularSave;
	Auxillary::alignedvector::ZV linDVscfRadial, linDVscfRadialSave;
	dataRegular_.resize(boost::extents[nR][numModes_][nrPC]);
	dataRadial_.resize(boost::extents[nR][numModes_][nAngChnlMax][nRadMax]);

	std::fill(dataRegular_.data(), dataRegular_.data()+ dataRegular_.size(), 0.0);
	for ( int atomIndex : primCell->get_atom_symmetry()->get_list_irreducible_atoms())
	{
		LatticeStructure::Symmetry const & siteSymmetry = primCell->get_site_symmetry(atomIndex);

		// obtain the pseudo inverse of the reducible set of displacements for this atom
		displCollection->generate_pseudo_inverse_reducible_displacements_for_atom(atomIndex, pseudoInverseDisplacements);

		// Symmetry-expand the set of potential variations belonging to this atom
		this->symmetry_expand_displacement_data(
				siteSymmetry,
				displCollection,
				atomIndex,
				nASC,
				potentialChange,
				deltaVRegularSymExpanded,
				deltaVRadialSymExpanded);
		assert(deltaVRegularSymExpanded.shape()[0] == displCollection->get_num_red_displacements_for_atom(atomIndex) );
		assert(deltaVRegularSymExpanded.shape()[1] == nrSC );
		assert(deltaVRadialSymExpanded.shape()[0] == displCollection->get_num_red_displacements_for_atom(atomIndex) );
		assert(deltaVRadialSymExpanded.shape()[1] == nASC );
		assert(deltaVRadialSymExpanded.shape()[2] == nAngChnlMax );
		assert(deltaVRadialSymExpanded.shape()[3] == nRadMax );

		// multiplying the Moore-Penrose pseudo inverse, of the displacements,
		// we obtain the least square fit of transpose of the linear displacement potential
		linDVscfRegular.resize(3*nrSC);
		linAlg.matrix_matrix_prod( pseudoInverseDisplacements, deltaVRegularSymExpanded, linDVscfRegular, 3, nrSC );
		// since radial data is complex, so must be the displacements for the call to BLAS
		pseudoInverseDisplacementsCmplx.resize(boost::extents[pseudoInverseDisplacements.shape()[0]][pseudoInverseDisplacements.shape()[1]]);
		std::copy_n(pseudoInverseDisplacements.data(), pseudoInverseDisplacements.size(),
				pseudoInverseDisplacementsCmplx.data());
		linDVscfRadial.resize(3 * nASC*nAngChnlMax*nRadMax);
		linAlg.matrix_matrix_prod( pseudoInverseDisplacementsCmplx, deltaVRadialSymExpanded, linDVscfRadial, 3, nASC*nAngChnlMax*nRadMax );

		// The lambda below has to match the layout of linDVscfRadial as defined by
		// the call to the matrix matrix product by which it was set.
		auto linRadialDataLayout = [&] (int iRadial, int iAngChnl, int scAtomIndex, int iDisplDirection) {
			assert((iAngChnl>=0)&&(iAngChnl < nAngChnlMax)
					&&(iDisplDirection>=0)&&(iDisplDirection < 3)
					&&(iRadial>=0)&&(iRadial < nRadMax)
					&&(scAtomIndex>=0)&&(scAtomIndex < nASC));
			// The layout is a 3 x (nASC,nAngChnlMax,nRadMax) matrix, or better tensor, in c-matrix-layout.
			// This is confirmed by the checks after the call to symmetry_expand_displacement_data by the multiarray
			// deltaVRadialSymExpanded. The matrix multiplication effectively replaces the displacement dimension by 3.
			return iRadial + nRadMax*(iAngChnl + nAngChnlMax*(scAtomIndex+nASC*iDisplDirection));
		};

		// Here we define a helper lambda that defines the layout as required by symOp.rotate_radial_vector_data
		// The layout as required by symOp.rotate_radial_vector_data is that 1) any size'd block of data is transformed
		// for 2) a given l,m channel number as the before slowest and 3) the spatial direction index as the slowest running
		// dimension. The requires copying around data which is done in the following.
		auto rotate_radial_vector_data_layout = [&] (int iRadial, int iAngChnl, int scAtomIndex, int iDisplDirection){
			return iRadial + nRadMax*(scAtomIndex+nASC*(iAngChnl + nAngChnlMax*iDisplDirection));
		};

		linDVscfRadialSave.assign(linDVscfRadial.begin(), linDVscfRadial.end());
		for (int iBlock = 0 ; iBlock < 3; ++iBlock)
			for (int iASC = 0 ; iASC < nASC; ++iASC)
				for (int iAngChnl = 0 ; iAngChnl < nAngChnlMax; ++iAngChnl)
					std::copy_n(&linDVscfRadialSave[linRadialDataLayout(0, iAngChnl, iASC, iBlock)], nRadMax,
							&linDVscfRadial[rotate_radial_vector_data_layout(0, iAngChnl, iASC, iBlock)]);

		// go through the star of this atom and set the internal data
		// note that the fitted data is now a vector, so it has to be transformed according
		// to the inverse symmetry operation
		for (std::pair<int,int> atomStar : primCell->get_atom_symmetry()->get_star_atom_indices(atomIndex))
		{
			const int atomIndexStar = atomStar.first;
			const int symIndex = atomStar.second;

			if ( symIndex != siteSymmetry.get_identity_index() )
			{
				// in case we are transformating the data, we make a copy for the next iteration
				linDVscfRegularSave.assign(linDVscfRegular.begin(), linDVscfRegular.end());
				linDVscfRadialSave.assign(linDVscfRadial.begin(), linDVscfRadial.end());
				auto symOp = siteSymmetry.get_sym_op(symIndex);
				symOp.transform_vector_field_regular_grid(supercellGrid, linDVscfRegular, /*transpose linDVscfRegular = */true);

				// Please note: here we rotate the data, but also the location, i.e. the atom
				//				the data refers to is affected by the transformation. This
				//				will be accounted for later, when we assign the data to the
				//				appropriate location.
				symOp.rotate_radial_vector_data(LMax, nRadMax*nASC, linDVscfRadial.begin(), linDVscfRadial.end());
			}

			// First get the regular part sorted.
			// Convert into the format of the linear derivative potential in real space, this also means
			// folding the data to correct primitive cell with the atom centered in the middle
			for ( int ir = 0 ; ir < nrSC; ++ir)
				for ( int xi_displ = 0 ; xi_displ < 3; ++xi_displ)
				{
					int modeIndex = Auxillary::memlayout::mode_layout(atomIndexStar,xi_displ);
					int ir1uc = perAtomGridIndexMap_[atomIndexStar][ir];
					int iR = this->R_vector_to_index( 	perAtomGridIndexMapRVector_[atomIndexStar][ir][0],
														perAtomGridIndexMapRVector_[atomIndexStar][ir][1],
														perAtomGridIndexMapRVector_[atomIndexStar][ir][2]);
					dataRegular_[iR][modeIndex][ir1uc] = linDVscfRegular[xi_displ*nrSC+ir];
				}

			// Now conclude with the radial part.
			for (int iASC = 0 ; iASC < nASC; ++iASC)
			{
				// Above, we have only applied the rotation to a given atomic site, but left the position
				// in the list of atoms unaffected.
				// This locates the atom after application of the symmetry operation.
				int iASCRot = superCell->get_atom_symmetry()->atom_rot_map(symIndex, iASC);
				int atomIndexPrimitive = primToSCCon->get_equiv_atom_primitive(iASCRot);

				std::array<int,3> latticeVectorPrimitiveToSC =
						primToSCCon->get_supercell_vector(primToSCCon->get_equiv_atom_primitive_lattice_vector_index(iASCRot));
				// possibly map this vector back to the embedding cell
				this->lattice_vector_map_back_embedding(latticeVectorPrimitiveToSC);
				int iR = this->R_vector_to_index(latticeVectorPrimitiveToSC[0], latticeVectorPrimitiveToSC[1], latticeVectorPrimitiveToSC[2]);
				for ( int xi_displ = 0 ; xi_displ < 3; ++xi_displ)
				{
					int modeIndex = Auxillary::memlayout::mode_layout(atomIndexPrimitive, xi_displ);
					for ( int ilm = 0 ; ilm < nAngChnlMax; ++ilm)
						for ( int iradial = 0 ; iradial < nRadMax; ++iradial)
							dataRadial_[iR][modeIndex][ilm][iradial] = linDVscfRadial[rotate_radial_vector_data_layout(iradial, ilm, iASC, xi_displ)];
				}
			}

			// swap back the copy
			if ( symIndex != siteSymmetry.get_identity_index() )
			{
				std::swap(linDVscfRegularSave, linDVscfRegular);
				std::swap(linDVscfRadialSave, linDVscfRadial);
			}
		}
	}

	//keep a copy for lattice information, symmetry ect ...
	primitiveCell_ = primCell;
	primitiveCellGrid_ = std::move(primcellGrid);

	this->clean_displacement_potential();
}

void
DisplacementPotential::compute_dvscf_q(
		std::vector<double> const & qVect,
		Auxillary::Multi_array<double,2> const & modes,
		Auxillary::Multi_array<std::complex<double>,3> const & dynamicalMatrices,
		std::vector<double> const & masses,
		Auxillary::alignedvector::CV & dvscf,
		std::vector<Auxillary::alignedvector::CV> & buffer,
		double freqCutoff) const
{
	assert( qVect.size() % 3 == 0 );
	const int nq = qVect.size()/3;
	const int nR = this->get_num_R();
	const int nAUC = numModes_/3;
	assert( modes.size() == numModes_*nq );
	assert( (dynamicalMatrices.size()/nq) / (numModes_*numModes_) == 1 );
	assert( (masses.size()*3) / numModes_ == 1 );
	freqCutoff = freqCutoff < 1e-5 ? 1e-5 : freqCutoff;

	// introduce names for the various buffers
	buffer.resize(2);
	Auxillary::alignedvector::CV & ftDisplPot = buffer[0];
	ftDisplPot.resize(numModes_*nptsRealSpace_);
	Auxillary::alignedvector::CV & dynmatMassBuffer = buffer[1];

	dvscf.assign(nq*numModes_*nptsRealSpace_, 0.0f);
	Algorithms::LinearAlgebraInterface linAlg;
	for ( int iq = 0 ; iq < nq; ++iq)
	{
		// this computes the Fourier transform. Note that internally
		// the data is stored with each atom in the center. After the Fourier
		// transform, we want to be back in the reference frame of the original
		// primitive cell. Thus we have to shift back by adding the atomic coordinate to each Lattice vector.
		std::fill(ftDisplPot.begin(), ftDisplPot.end(), std::complex<float>(0));
		for (int ia = 0 ; ia < nAUC; ++ia)
		{
			double tau[3] = {primitiveCell_->get_atoms_list()[ia].get_position()[0],
					primitiveCell_->get_atoms_list()[ia].get_position()[1],
					primitiveCell_->get_atoms_list()[ia].get_position()[2] };
			for ( int iR = 0 ; iR < nR; ++iR)
			{
				float dprod = -2.0*M_PI*(qVect[iq*3+0]*(RVectors_[iR][0] + tau[0])
										+qVect[iq*3+1]*(RVectors_[iR][1] + tau[1])
										+qVect[iq*3+2]*(RVectors_[iR][2] + tau[2]));
				std::complex<float> phase = std::exp( std::complex<float>(0,dprod));

				for (int iDispl = 0 ; iDispl < 3; ++iDispl)
				{
					int mu = Auxillary::memlayout::mode_layout(ia,iDispl);
					for (int ir = 0 ; ir < nptsRealSpace_; ++ir)
						ftDisplPot[mu*nptsRealSpace_+ir] += phase*dataRegular_[iR][mu][ir];
				}
			}
		}

		// Here we multiply the dynamical matrix to the Fourier transformed data
		dynmatMassBuffer.resize(numModes_*numModes_);
		for ( int mu1 = 0 ; mu1 < numModes_; ++mu1 )
			for ( int mu2 = 0 ; mu2 < numModes_; ++mu2 )
				dynmatMassBuffer[mu1*numModes_+mu2] =
						std::complex<float>(dynamicalMatrices[iq][mu1][mu2]) /  std::sqrt(float(masses[mu1/3]));

		auto r_ptr = &dvscf[iq*numModes_*nptsRealSpace_];
		linAlg.call_gemm( 't', 'n',
				numModes_, nptsRealSpace_ , numModes_,
				std::complex<float>(1.0f), dynmatMassBuffer.data(), numModes_,
				ftDisplPot.data(), nptsRealSpace_,
				std::complex<float>(0.0f), r_ptr, nptsRealSpace_);
	}

	// multiply the phase exp(i q*r) to make it lattice periodic
	this->multiply_phase_qr(qVect, dvscf);

	// symmetrize the lattice periodic displacement potential
//	this->symmetrize_periodic_dvscf_q(qVect, dvscf);

	// go to units of eV by dividing by the frequency and multiplying the conversion factor
	for ( int iq = 0 ; iq < nq; ++iq)
		for ( int mu = 0 ; mu < numModes_; ++mu )
		{
			float scale = (modes[iq][mu] < freqCutoff ? 0.0f :
					static_cast<float>(Auxillary::units::SQRT_HBAR_BY_2M_THZ_TO_ANGSTROEM / std::sqrt(modes[iq][mu])) );
			for (int ir = 0 ; ir < nptsRealSpace_; ++ir)
			{
				dvscf[(iq*numModes_+mu)*nptsRealSpace_+ir] *= scale;
				assert(dvscf[(iq*numModes_+mu)*nptsRealSpace_+ir] == dvscf[(iq*numModes_+mu)*nptsRealSpace_+ir]);
			}
		}
}

void
DisplacementPotential::query_q(
		std::vector<double> const & qVect,
		std::map<LatticeStructure::RegularBareGrid::GridPoint, std::vector<int>> & qGridCorseGainToQIndex) const
{
	assert(qVect.size()%3 == 0);
	qGridCorseGainToQIndex.clear();
	const int nq = qVect.size()/3;
	const double gPrec = primitiveCell_->get_symmetry().get_symmetry_prec();

	for (int iq = 0 ; iq < nq ;++iq)
	{
		LatticeStructure::RegularBareGrid::GridPoint g(std::move(std::vector<double>(&qVect[iq*3], &qVect[iq*3]+3)), gPrec);
		auto ret = qGridCorseGainToQIndex.insert(std::make_pair(std::move(g), std::vector<int>()));
		ret.first->second.push_back(iq);
	}
}

int
DisplacementPotential::get_num_R() const
{
	return RVectors_.shape()[0];
}

int
DisplacementPotential::get_num_modes() const
{
	return numModes_;
}

int
DisplacementPotential::mem_layout(int ir, int mu, int iR ) const
{
	assert(ir < this->nptsRealSpace_ );
	assert(mu < this->numModes_);
	assert(iR < this->get_num_R());
	return ir+nptsRealSpace_*(mu+numModes_*iR);
}

LatticeStructure::RegularBareGrid const &
DisplacementPotential::get_real_space_grid() const
{
	return primitiveCellGrid_;
}

void
DisplacementPotential::write_dvscf(
		int atomIndex, int xi,
		std::string filename) const
{
	assert((atomIndex >= 0) && (atomIndex<numModes_/3));
	assert((xi >= 0) && (xi < 3));
	assert(this->get_num_R() > 0);
	assert(perAtomGridIndexMap_.shape()[0]>atomIndex);

	const int mu = Auxillary::memlayout::mode_layout(atomIndex, xi);

	IOMethods::WriteVASPRealSpaceData writer;
	std::string comment = "Displacement potential in real space; atom # "+ std::to_string(atomIndex) +", vibration in dir. "
			+ (xi == 0 ? "x" : (xi == 1 ? "y" : "z") ) + "\n";

	LatticeStructure::UnitCell supercell;

	const int nASC = primitiveCell_->get_atoms_list().size()*RVectors_.shape()[0];

	// build a supercell where the displaced atom is in the center.
	auto atomsUC = primitiveCell_->get_atoms_list();
	auto center = atomsUC[atomIndex].get_position();
	primToSuperCell_->primitive_to_supercell_coordinates(center);

	std::vector<LatticeStructure::Atom> scAtoms;
	scAtoms.reserve(nASC);
	for (int iR = 0 ; iR < RVectors_.shape()[0]; ++iR)
	{
		for (auto a: atomsUC)
		{
			auto pos = a.get_position();
			for (int i = 0 ; i < pos.size(); ++i)
				pos[i] = (RVectors_[iR][i]+pos[i]);
			primToSuperCell_->primitive_to_supercell_coordinates(pos);
			for (int i = 0 ; i < pos.size(); ++i)
				pos[i] = pos[i] - center[i];
			a.set_position(std::move(pos));
			scAtoms.push_back(std::move(a));
		}
	}
	supercell.initialize(	std::move(scAtoms),
							primitiveCell_->get_lattice().build_supercell(primToSuperCell_->get_supercell_matrix()),
							primitiveCell_->get_symmetry());

	std::vector<double> dataSC(perAtomGridIndexMap_.shape()[1]);
	for ( int irSC = 0 ; irSC < perAtomGridIndexMap_.shape()[1]; ++irSC )
	{
		int iR = this->R_vector_to_index(	perAtomGridIndexMapRVector_[atomIndex][irSC][0],
											perAtomGridIndexMapRVector_[atomIndex][irSC][1],
											perAtomGridIndexMapRVector_[atomIndex][irSC][2]);
		dataSC[irSC] = dataRegular_[iR][mu][perAtomGridIndexMap_[atomIndex][irSC]];
	}

	writer.write_file(	filename,
						comment,
						embeddingSuperCellGrid_.get_grid_dim(),
						std::make_shared<decltype(supercell)>(std::move(supercell)),
						dataSC );
}

void
DisplacementPotential::write_dvscf_q(
		std::vector<double> const & qVect,
		std::vector<int> modeIndices,
		Auxillary::Multi_array<double,2> const & modes,
		Auxillary::Multi_array<std::complex<double>,3> const & dynamicalMatrices,
		std::vector<double> const & masses,
		std::string filename) const
{
	assert( qVect.size()%3 == 0 );
	int nq = qVect.size()/3;

	Auxillary::alignedvector::CV qDataPlus;
	std::vector<Auxillary::alignedvector::CV> buffers;
	this->compute_dvscf_q(qVect, modes, dynamicalMatrices, masses, qDataPlus, buffers);

	int nr = primitiveCellGrid_.get_num_points();
	assert(qDataPlus.size() == nq*numModes_*nr);

	auto mod_filename = [] (std::string filename, int iq, int mu) {
		if ( filename.length() >= 5 )
			if ( filename.substr(filename.length()-4,filename.length()).compare( ".dat" ) == 0 )
				filename.insert(filename.length()-4,"_"+std::to_string(iq)+"_"+std::to_string(mu));
		return filename;
	};

	IOMethods::WriteVASPRealSpaceData writer;
	for ( int iq = 0; iq < nq ; ++iq )
	{
		for ( auto mu : modeIndices )
		{
			assert( (mu >= 0) && (mu < numModes_) );
			std::string comment = "Lattice periodic part of the displacement potential dvscf(r)/du(q,mu); q = ("
					+std::to_string(qVect[iq*3+0])+" , "+std::to_string(qVect[iq*3+1])+" , "+std::to_string(qVect[iq*3+2])
					+") mode # "+ std::to_string(mu) +". Units are eV.\n";
			std::vector<double> data(nr);
			for ( int ir = 0 ; ir < nr ; ++ir )
				data[ir] = std::real(qDataPlus[(iq*numModes_+mu)*nr+ir]);
			writer.write_file( mod_filename(filename,iq,mu), comment,
					primitiveCellGrid_.get_grid_dim(), primitiveCell_, data );
		}
	}
}

void
DisplacementPotential::clean_displacement_potential()
{
	const int numAtoms = this->get_num_modes()/3;
	const int nr = primitiveCellGrid_.get_num_points();
	const int nR = this->get_num_R();

	// apply the sum rule
	for ( int xi = 0 ; xi < 3 ; ++xi)
	{
		double integral = 0.0;
		for (int iR = 0; iR < nR ; ++iR)
			for (int ia = 0 ; ia < numAtoms; ++ia)
				for (int ir = 0; ir < nr ; ++ir)
					integral += dataRegular_[iR][Auxillary::memlayout::mode_layout(ia,xi)][ir];
		integral /= static_cast<double>(nr*nR*numAtoms);

		for (int iR = 0; iR < nR ; ++iR)
			for (int ia = 0 ; ia < numAtoms; ++ia)
				for (int ir = 0; ir < nr ; ++ir)
					dataRegular_[iR][Auxillary::memlayout::mode_layout(ia,xi)][ir] -= integral;
	}
}

void
DisplacementPotential::symmetry_expand_displacement_data(
		LatticeStructure::Symmetry const & siteSymmetry,
		std::shared_ptr<const LatticeStructure::AtomDisplacementCollection> displColl,
		int atomIndex,
		int numAtomsSupercell,
		std::vector<std::shared_ptr<const PotentialChangeIrredDisplacement>> const & potentialChange,
		Auxillary::Multi_array<double,2> & deltaVRegularSymExpanded,
		Auxillary::Multi_array<std::complex<double>,4> & deltaVRadialSymExpanded) const
{
	const int nRD = displColl->get_num_red_displacements_for_atom(atomIndex);
	const int nrSC = std::distance(potentialChange[0]->begin_regular_data(), potentialChange[0]->end_regular_data());
	deltaVRegularSymExpanded.resize(boost::extents[nRD][nrSC]);

	const int nRadMax = potentialChange[0]->get_max_num_radial_elements();
	const int nAngChnlMax = potentialChange[0]->get_max_num_angular_moment_channels();
	deltaVRadialSymExpanded.resize(boost::extents[nRD][numAtomsSupercell][nAngChnlMax][nRadMax]);

	auto starDisplacements = displColl->get_symmetry_relation_red_displacements_for_atom(atomIndex);
	for (int redDispl = 0 ; redDispl < starDisplacements.size() ; ++redDispl  )
	{
		const int irredDisplIndex = starDisplacements[redDispl].first;
		const int irredToRedSymIndex = starDisplacements[redDispl].second;
		assert((irredDisplIndex >= 0) && (irredDisplIndex < potentialChange.size()));

		auto transformPotentialChange = potentialChange[irredDisplIndex];
		if (irredToRedSymIndex != siteSymmetry.get_identity_index())
		{// for non idenity symmetry replace the ptr with a transformed one ...
			PotentialChangeIrredDisplacement potentialChangeCopy = *potentialChange[irredDisplIndex];
			potentialChangeCopy.transform(siteSymmetry.get_sym_op(irredToRedSymIndex));
			transformPotentialChange = std::make_shared<const PotentialChangeIrredDisplacement>(std::move(potentialChangeCopy));
		}

		assert(std::distance(transformPotentialChange->begin_regular_data(),transformPotentialChange->end_regular_data())
			== nrSC);
		std::copy(transformPotentialChange->begin_regular_data(),
				transformPotentialChange->end_regular_data(),
				&deltaVRegularSymExpanded[redDispl][0]);
		for (int iaSC = 0 ; iaSC < numAtomsSupercell; ++iaSC)
		{
			assert(std::distance(transformPotentialChange->begin_radial_data(iaSC),transformPotentialChange->end_radial_data(iaSC))
				== nAngChnlMax*nRadMax);
			std::copy(transformPotentialChange->begin_radial_data(iaSC),
					transformPotentialChange->end_radial_data(iaSC),
					&deltaVRadialSymExpanded[redDispl][iaSC][0][0]);
		}
	}
}

void
DisplacementPotential::symmetrize_periodic_dvscf_q(
		std::vector<double> const & qpoints,
		Auxillary::alignedvector::CV & dvscfData) const
{
	const int nq = qpoints.size()/3;
	const int nM = this->get_num_modes();
	assert(dvscfData.size()==nq*nM*nptsRealSpace_);

	LatticeStructure::Symmetry const & fullSymmetryGroup = primitiveCell_->get_symmetry();
	std::vector<std::vector<int>> rotMap;
	Algorithms::compute_grid_rotation_map_no_shift(
			primitiveCellGrid_,
			fullSymmetryGroup,
			rotMap);

	Auxillary::alignedvector::CV rotBuffer(3*3);
	for ( int iq = 0 ; iq < nq ; ++iq)
	{
		std::vector<double> q(&qpoints[iq*3], &qpoints[iq*3]+3);
		auto smallGroupQ = fullSymmetryGroup;
		smallGroupQ.set_reciprocal_space_sym(true);
		auto dropedSymmetries = smallGroupQ.small_group(q);
		const int nsymQ = smallGroupQ.get_num_symmetries_no_T_rev();

		// for the pure identity symmetry there is nothing to be done.
		if (nsymQ == 1)
			continue;

		std::set<int> dropSym(dropedSymmetries.begin(), dropedSymmetries.end());

		for (int inu = 0 ; inu < nM; ++inu)
		{
			Auxillary::alignedvector::CV dvscfToSym(&dvscfData[(iq*nM+inu)*nptsRealSpace_],
													&dvscfData[(iq*nM+inu)*nptsRealSpace_]+nptsRealSpace_);
			for (int isym = 0 ; isym < fullSymmetryGroup.get_num_symmetries(); ++isym)
			{
				if ( dropSym.find(isym) != dropSym.end())
					continue;

				// Identity symmetry is included with the copy above.
				if (isym == fullSymmetryGroup.get_identity_index())
					continue;

				for (int ir = 0 ; ir < nptsRealSpace_; ++ir)
				{
					const int irRot = rotMap[isym][ir];
					dvscfToSym[irRot] += dvscfData[(iq*nM+inu)*nptsRealSpace_+ir];
				}
			}
			for (int ir = 0 ; ir < nptsRealSpace_; ++ir)
			{
				dvscfData[(iq*nM+inu)*nptsRealSpace_+ir] = dvscfToSym[ir]/static_cast<float>(nsymQ);
			}
		}
	}
}

void
DisplacementPotential::multiply_phase_qr(
		std::vector<double> const & qpoints,
		Auxillary::alignedvector::CV & dvscf_q) const
{
	const int nq = qpoints.size()/3;
	const int nM = this->get_num_modes();
	assert(dvscf_q.size()==nq*nM*nptsRealSpace_);

	// the code below attempts to reduce expensive complex phase calculations in this performance critical part of the code
	const int nrx = primitiveCellGrid_.get_grid_dim()[0];
	const int nry = primitiveCellGrid_.get_grid_dim()[1];
	const int nrz = primitiveCellGrid_.get_grid_dim()[2];

	for (int iq = 0 ; iq < nq ; ++iq)
	{
		auto phaseIncX = std::complex<double>( std::cos(2.0*M_PI*qpoints[iq*3+0]/static_cast<double>(nrx)),
											   std::sin(2.0*M_PI*qpoints[iq*3+0]/static_cast<double>(nrx)) );
		auto phaseIncY = std::complex<double>( std::cos(2.0*M_PI*qpoints[iq*3+1]/static_cast<double>(nry)),
											   std::sin(2.0*M_PI*qpoints[iq*3+1]/static_cast<double>(nry)) );
		auto phaseIncZ = std::complex<double>( std::cos(2.0*M_PI*qpoints[iq*3+2]/static_cast<double>(nrz)),
											   std::sin(2.0*M_PI*qpoints[iq*3+2]/static_cast<double>(nrz)) );

		#ifndef NDEBUG
		auto rVec = primitiveCellGrid_.get_all_vectors_grid();
		#endif

		for ( int mu = 0 ; mu < numModes_; ++mu )
		{
			std::complex<double> totalPhaseGP(1.0);
			for (int iz = 0 ; iz < nrz; ++iz)
			{
				if (iz == nrz/2)
					totalPhaseGP *= std::complex<double>(  std::cos(2.0*M_PI*qpoints[iq*3+2]),
														  -std::sin(2.0*M_PI*qpoints[iq*3+2]) );
				for (int iy = 0 ; iy < nry; ++iy)
				{
					if (iy == nry/2)
						totalPhaseGP *= std::complex<double>(  std::cos(2.0*M_PI*qpoints[iq*3+1]),
															  -std::sin(2.0*M_PI*qpoints[iq*3+1]) );
					for (int ix = 0 ; ix < nrx; ++ix)
					{
						if (ix == nrx/2)
							totalPhaseGP *= std::complex<double>(  std::cos(2.0*M_PI*qpoints[iq*3+0]),
																  -std::sin(2.0*M_PI*qpoints[iq*3+0]) );
						const int ir = ix+nrx*(iy+nry*iz);
#ifndef NDEBUG
						// in debug mode, we verify that the phase (computed directly) equal the stepping
						// algorithm used here to avoid expensive trigonometric functions.
						float arg = 2.0*M_PI*(  qpoints[iq*3+0]*rVec[ir*3+0]
											  + qpoints[iq*3+1]*rVec[ir*3+1]
											  +	qpoints[iq*3+2]*rVec[ir*3+2]);
						auto p = std::complex<float>(std::cos(arg), std::sin(arg));
						assert(std::abs(std::complex<float>(totalPhaseGP)-p)<1e-6);
#endif
						dvscf_q[(iq*numModes_+mu)*nptsRealSpace_+ir] *= static_cast<std::complex<float>>(totalPhaseGP);
						totalPhaseGP *= phaseIncX;
					}
					totalPhaseGP *= phaseIncY;
				}
				totalPhaseGP *= phaseIncZ;
			}
		}
	}
}

int
DisplacementPotential::R_vector_to_index(int Rx, int Ry, int Rz) const
{
	assert((Rx>=-embeddingRVectorBox_[0]) && (Rx<=embeddingRVectorBox_[0]));
	assert((Ry>=-embeddingRVectorBox_[1]) && (Ry<=embeddingRVectorBox_[1]));
	assert((Rz>=-embeddingRVectorBox_[2]) && (Rz<=embeddingRVectorBox_[2]));
	int iR = RVectoToIndexMap_[Rx][Ry][Rz];
	assert((iR >= 0) && (iR < this->get_num_R()));
	return iR;
}

void
DisplacementPotential::lattice_vector_map_back_embedding(std::array<int,3> & latticeVector) const
{
	for (int i = 0 ; i < 3 ; ++i)
	{
		if ( latticeVector[i] < -embeddingRVectorBox_[i] )
			latticeVector[i] += (2*embeddingRVectorBox_[i]+1);
		if ( latticeVector[i] > embeddingRVectorBox_[i] )
			latticeVector[i] -= (2*embeddingRVectorBox_[i]+1);
		assert((latticeVector[i]>=-embeddingRVectorBox_[i]) && (latticeVector[i]<=embeddingRVectorBox_[i]));
	}
}

Auxillary::Multi_array<float,3> const &
DisplacementPotential::get_data_regular_grid() const
{
	return dataRegular_;
}

Auxillary::Multi_array<std::complex<float>,4> const &
DisplacementPotential::get_data_radial_grid() const
{
	return dataRadial_;
}

} /* namespace PhononStructure */
} /* namespace elephon */
