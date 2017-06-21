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
#include <assert.h>
#include <cmath>
#include <set>
#include <map>
#include <iostream>
#include "boost/multi_array.hpp"

namespace elephon
{
namespace PhononStructure
{

void
ForceConstantMatrix::build(  LatticeStructure::UnitCell unitCell,
		LatticeStructure::UnitCell const & superCell,
		std::vector<LatticeStructure::AtomDisplacement> const & irredDispl,
		std::vector<std::vector<double>> forces)
{
	//Forces F(i,xi) acting on atom i in direction xi are given by
	//	F(i,xi) = \Phi(i,j,xi,xj) * u(j,xi)
	//where \Phi is the matrix of force constants and u are the displacements of an atom j in direction xj.
	//Here we solve for given (input) F and u for Phi.
	//
	//The general approach is the following (compare https://doi.org/10.1103/PhysRevLett.78.4063)
	//
	//	1)	Construct the force constants at each irreducible displacement by writing
	//		 F(i,xi) = U(k) * \Phi(k,i,xi)
	this->set_supercell_dim( unitCell, superCell );
	numModes_ = unitCell.get_atoms_list().size()*3;

	assert( forces.size() == (irredDispl.size() + 1) );

	//Regenerate the list of irreducible atoms. This is needed to map the input forces
	// which come for non-equivalent displacements at irreducible atoms into the reducible set.
	std::vector<LatticeStructure::Atom> irredAtoms;
	std::vector<int> redToIrredAtoms, symRedToIrredAtoms;
	std::vector< std::vector<int> > irredToRedAtoms, symIrredToRedAtoms;
	LatticeStructure::SymmetryReduction<LatticeStructure::Atom>(
			unitCell.get_symmetry(),
			unitCell.get_atoms_list(),  irredAtoms,
			redToIrredAtoms, symRedToIrredAtoms,
			irredToRedAtoms, symIrredToRedAtoms);

	//enforce sum F=0
	this->drift_clean_forces( forces );

	//define a consistent atom numbering in the primitive cell.
	//We choose the unitCell.get_atoms_list() order
	std::map< LatticeStructure::Atom, int > primitiveCellLookup;
	for ( int ia = 0 ; ia < unitCell.get_atoms_list().size(); ++ia )
		primitiveCellLookup.insert( std::move( std::make_pair( unitCell.get_atoms_list()[ia], ia ) ) );

	int naSC = superCell.get_atoms_list().size();

	Algorithms::LinearAlgebraInterface linAlg;

	//For each such atom, solve from the reducible displacements for the
	//local matrix of force constants
	std::vector<double> irreducibleMatrixOfForceConstants(
			irredAtoms.size()*superCell.get_atoms_list().size()*9);
	int iForceField = 0;
	for ( int ia = 0 ; ia < irredAtoms.size(); ++ia )
	{
		//regenerate displacements to make the connection with the input forces
		LatticeStructure::Symmetry const & sSym = unitCell.get_site_symmetry(ia);

		std::vector<LatticeStructure::AtomDisplacement> irredNotUsed,reducibleNotUsed;
		bool symmetricDispl = not irredDispl[iForceField].is_plus_minus_displ_equivalent();
		std::vector<int> redToIrredDispl, symRedToIrredDispl;
		std::vector< std::vector<int> > irredToRedDispl, symIrredToRedDispl;
		unitCell.get_site_displacements( irredAtoms[ia],
				symmetricDispl, sSym,
				1, // not used here
				irredNotUsed,reducibleNotUsed,
				redToIrredDispl, symRedToIrredDispl,
				irredToRedDispl, symIrredToRedDispl);

		//Otherwise the system is underdetermined - this is not supposed to happen.
		int nRedDispl = redToIrredDispl.size();
		assert( nRedDispl >= 3 );

		//build the transpose of the U matrix from the irreducible displacements
		std::vector<double> U( 3 * nRedDispl );
		for ( int idir = 0 ; idir < irredToRedDispl.size();  ++idir)
		{
			auto d = irredDispl[iForceField+idir].get_direction();
			for ( auto &xi : d )
				xi *= irredDispl[iForceField].get_magnitude();
			for ( int istar = 0 ; istar < symIrredToRedDispl[idir].size(); ++istar)
			{
				auto drot = d;
				sSym.rotate_cartesian( symIrredToRedDispl[idir][istar], drot.begin(), drot.end() );
				std::copy( drot.begin(), drot.end(),  &U[irredToRedDispl[idir][istar]*3] );
			}
		}

		//compute the pseudo inverse
		std::vector<double> piU;
		linAlg.pseudo_inverse( std::move(U), nRedDispl, 3, piU );

		//Obtain the displaced atom in the coordinates of the supercell
		LatticeStructure::Atom displacedAtomInSupercell = irredAtoms[ia];
		auto p = irredAtoms[ia].get_position();
		displacedAtomInSupercell.set_position( {p[0]/double(supercellDim_[0]),
												p[1]/double(supercellDim_[1]),
												p[2]/double(supercellDim_[2]) } );

		//For reducible displacements in the star of an irreducible one, this means, we need to
		//	1) identify position g_\alpha^-1 * i, where i is a supercell atom index and g_\alpha is a rotation
		//		around the current site index.
		//	2) rotate the 3-vector with the forces by g_\alpha
		//For all reducible displacement at this site, build the transpose of the matrix F
		std::map<LatticeStructure::Atom,int> rotLookup;
		this->shift_atoms_into_map(
				superCell.get_atoms_list(), displacedAtomInSupercell,
				rotLookup);
		//rotMapsSiteSymmetry tells for site symmetry where each atom is
		//transported while going from the irreducible to the reducible displacement
		std::vector< std::vector<int> > rotMapsSiteSymmetry( nRedDispl );
		for ( int id = 0 ; id < nRedDispl; ++id)
		{
			int isym = symRedToIrredDispl[id];
			if ( isym != sSym.get_identity_index() )
			{
				int isymInv = sSym.get_index_inverse( isym ); //because we go from irreducible to reducible
				this->transform_map(rotLookup, sSym.get_sym_op(isym), rotMapsSiteSymmetry[id]);
			}
			else
			{
				rotMapsSiteSymmetry[id].resize( superCell.get_atoms_list().size() );
				for ( int i = 0; i < rotMapsSiteSymmetry[id].size(); ++i )
					rotMapsSiteSymmetry[id][i] = i;
			}
		}

		std::vector<double> matrixForceConstants(3*3);
		std::vector<double> F( 3*nRedDispl );
		for ( int iaS = 0 ; iaS < naSC ; ++iaS )
		{
			for ( int id = 0 ; id < nRedDispl; ++id)
			{
				int irredForceField = redToIrredDispl[id];
				if( iForceField+irredForceField >= forces.size() )
					throw std::runtime_error("Set of forces and number of irreducible atoms does not match");

				int iForceAtom = rotMapsSiteSymmetry[id][iaS];
				if( ! ( (iForceAtom >=0) && (iForceAtom*3 < forces[iForceField+irredForceField].size())) )
					throw std::runtime_error("Size of the force field data insufficient");

				std::copy( &forces[iForceField+irredForceField][iForceAtom*3],
						&forces[iForceField+irredForceField][iForceAtom*3]+3,
						&F[id*3] );

				//note, we have to use the inverse index, because we are going from reducible to irreducible
				sSym.rotate_cartesian( sSym.get_index_inverse(symRedToIrredDispl[id]),
						F.begin()+id*3, F.begin()+(id+1)*3 );
			}

			//multiplying the Moore-Penrose pseudo inverse, of the displacements,
			//we obtain the least square fit of transpose of the force constant matrix
			linAlg.matrix_matrix_prod( piU, F, matrixForceConstants, 3, 3 );

			//transpose back
			auto localMat = irreducibleMatrixOfForceConstants.data()+(ia*naSC+iaS)*9;
			for ( int xi = 0 ; xi < 3; ++xi)
				for ( int xj = 0 ; xj < 3; ++xj)
					localMat[xi*3+xj] = matrixForceConstants[xj*3+xi];
		}
		iForceField += irredToRedDispl.size();
	}

	//Create a table which maps a supercell atom index into a primitive cell coordinate
	//and a lattice vector
	std::vector< std::pair<int,std::vector<int> > > atomSuperCellToPrimitve( naSC );
	for ( int iaS = 0 ; iaS < naSC ; ++iaS )
	{
		//Locate the place where to store the data in the full array
		//find the position of this atom in units of the primitive cell
		auto atom = superCell.get_atoms_list()[iaS];
		std::vector<double> posAtomUnitsPrimitve = atom.get_position();
		std::vector<int> R(3);
		std::vector<double> tau(3);
		for ( int i = 0 ; i < 3 ; ++i )
		{
			//scale to units of the primitive cell
			posAtomUnitsPrimitve[i] *= supercellDim_[i];
			//We add a tiny bit to the vector so that 1.9999999 will not (incorrectly) map
			//to R = 1 and then 0.9999999 will (correctly) map to the lattice site at 0.
			R[i] = std::floor(posAtomUnitsPrimitve[i]+0.5+2*unitCell.get_symmetry().get_symmetry_prec());
			tau[i] = posAtomUnitsPrimitve[i]-R[i];
		}
		atom.set_position(tau);
		auto it = primitiveCellLookup.find( atom );
		if ( it == primitiveCellLookup.end() )
			throw std::runtime_error( "Problem locating supercell atom after mapping back to the primitive cell" );
		atomSuperCellToPrimitve[iaS] = std::move(std::make_pair( it->second, std::move(R) ) );
	}

	//Fill the reducible atomic sites using the equation
	//	C( beta1 T, beta2 0 ) = g . C( S(alpha1,T), S(alpha2,0) ) . g^-1
	//where S is a symmetry operation connecting alpha(1,2) with beta(1,2) and g is the point group part.

	std::map< LatticeStructure::Atom, int > superCellLookup;
	for ( int ia = 0 ; ia < superCell.get_atoms_list().size(); ++ia )
		superCellLookup.insert(std::move( std::make_pair( superCell.get_atoms_list()[ia], ia ) ) );

	data_.resize( this->get_num_R()*std::pow(3*unitCell.get_atoms_list().size(),2) );
	RSupercellMultiplicity_.resize(this->get_num_R());
	std::vector<double> matForceSlice( 9*naSC );
	std::vector<int> rotMap;
	for ( int irA = 0; irA < irredToRedAtoms.size(); ++irA )
	{
		//loop the star of the irreducible atom
		for ( int istar = 0; istar < symIrredToRedAtoms[irA].size(); ++istar )
		{
			//transform the prelimnary matrix of force constants according
			//to the present symmetry operation which takes us from the irreducible
			//to the reducible atom
			int isym = symIrredToRedAtoms[irA][istar];
			if ( isym == unitCell.get_symmetry().get_identity_index() )
			{
				std::copy(&irreducibleMatrixOfForceConstants[irA*naSC*9],
						&irreducibleMatrixOfForceConstants[irA*naSC*9]+naSC*9,
						matForceSlice.begin() );
			}
			else
			{
				//Re-shuffle the matrices according to the symmetry operation beta = S( alpha )
				this->transform_map( superCellLookup, unitCell.get_symmetry().get_sym_op(isym), rotMap );
				for ( int iaSC = 0 ; iaSC < naSC ; ++iaSC )
				{
					int iredA1 = rotMap[iaSC];
					std::copy(&irreducibleMatrixOfForceConstants[(irA*naSC+iaSC)*9],
							  &irreducibleMatrixOfForceConstants[(irA*naSC+iaSC)*9]+9,
							&matForceSlice[iredA1*9] );
				}

				//transform the set of 3x3 matrices for the set of atoms
				unitCell.get_symmetry().rotate_matrix_cartesian( isym,
						matForceSlice.begin(),
						matForceSlice.end());
			}
			//Convert into the format of the dynamical matrix in real space
			for ( int iaSC = 0 ; iaSC < naSC ; ++iaSC )
			{
				int ia1 = atomSuperCellToPrimitve[iaSC].first;
				int iaDispl = irredToRedAtoms[irA][istar];
				std::vector<int> R = atomSuperCellToPrimitve[iaSC].second;
				//R is the vector in C(mu1 R, mu2) and measured from [0,supercellDim_[
				//The dynamical matrix for the inifinite lattice is with the atoms
				//in the atom where the force acts in the middle and symmetric cells up
				//to a given range included. This range is thus [-supercellDim_/2, supercellDim_/2]
				//for even supercellDim_ and [-(supercellDim_-1)/2, (supercellDim_+1)/2] for
				//odd supercellDim_
				for ( int i = 0 ; i < 3 ; ++i)
				{
					//take R to -R (plus periodic repetition)
					R[i] = (supercellDim_[i] - R[i])%supercellDim_[i];
					//Map R to the range [-supercellDim_/2,supercellDim_[i]/2[
					//or [-(supercellDim_-1)/2, (supercellDim_+1)/2[ for odd
					R[i] = R[i] < (supercellDim_[i]+supercellDim_[i]%2)/2 ? R[i] : R[i]-supercellDim_[i];
				}

				//If R is -supercellDim_/2 (or the equivalent for odd supercellDim_)
				//we need to duplicate the data to the other side of the range
				std::vector< std::vector<int> > allR(1);
				allR[0] = std::move(R);
				for ( int i = 0 ; i < 3 ; ++i)
					if ( allR[0][i] == -(supercellDim_[i]-supercellDim_[i]%2)/2 )
					{
						auto allRCpy = allR;
						for ( auto &R : allRCpy )
							R[i] = (supercellDim_[i]+supercellDim_[i]%2)/2;
						allR.insert(std::end(allR),std::begin(allRCpy),std::end(allRCpy));
					}

				for ( auto R : allR )
				{
					for ( int i = 0 ; i < 3 ; ++i)
						R[i] = R[i] < 0 ? R[i]+RVectorDim_[i]:R[i];
					int ir = this->RVectorLayout(R[2],R[1],R[0]);
					RSupercellMultiplicity_[ir] = allR.size();
					for ( int xj = 0 ; xj < 3; ++xj)
						for ( int xi = 0 ; xi < 3; ++xi)
							data_[this->mem_layout(ir,ia1*3+xj,iaDispl*3+xi)] = matForceSlice[iaSC*9+xj*3+xi];
				}
			}
		}
	}

	double check = 0;
	for ( auto m : RSupercellMultiplicity_ )
		check += 1.0/double(m);
	if ( std::abs(supercellDim_[0]*supercellDim_[1]*supercellDim_[2]-check) > 1e-8 )
		throw std::logic_error("R grid weights do not add up to the supercell volume");

	tau_.resize( 3*unitCell.get_atoms_list().size() );
	for ( int ia = 0; ia < unitCell.get_atoms_list().size(); ++ia )
	{
		auto pos = unitCell.get_atoms_list()[ia].get_position();
		std::copy(pos.begin(), pos.end(), &tau_[3*ia]);
	}

	//Symmetrize
//	//enforce that this is a symmetric matrix
//	for ( int iRz = 0 ; iRz < RVectorDim_[2]; ++iRz )
//		for ( int iRy = 0 ; iRy < RVectorDim_[1]; ++iRy )
//			for ( int iRx = 0 ; iRx < RVectorDim_[0]; ++iRx )
//				for ( int mu1 = 0 ; mu1 < numModes_; ++mu1)
//					for ( int mu2 = mu1+1 ; mu2 < numModes_; ++mu2)
//					{
//						int miRx = RVectorDim_[0]-iRx-1;
//						int miRy = RVectorDim_[1]-iRy-1;
//						int miRz = RVectorDim_[2]-iRz-1;
//						double tmp = 0.5*(data_[ this->mem_layout(iRz,iRy,iRx,mu2,mu1)]
//								+data_[ this->mem_layout(miRz,miRy,miRx,mu1,mu2)] );
//						data_[ this->mem_layout(iRz,iRy,iRx,mu2,mu1) ] = tmp;
//						data_[ this->mem_layout(miRz,miRy,miRx,mu1,mu2)] = tmp;
//					}
}

double
ForceConstantMatrix::operator() (int Rz, int Ry, int Rx, int mu1, int mu2) const
{
	return data_[ this->mem_layout(Rz,Ry,Rx,mu2,mu1) ];
}

int
ForceConstantMatrix::RVectorLayout(int iRz, int iRy, int iRx ) const
{
	assert( (iRx >= 0) && (iRx < RVectorDim_[0]) );
	assert( (iRy >= 0) && (iRy < RVectorDim_[1]) );
	assert( (iRz >= 0) && (iRz < RVectorDim_[2]) );
	return (iRz*RVectorDim_[1]+iRy)*RVectorDim_[0]+iRx;
}

int
ForceConstantMatrix::mem_layout(int Rz, int Ry, int Rx, int mu1, int mu2) const
{
	return this->mem_layout(this->RVectorLayout(Rz,Ry,Rx),mu2,mu1);
};


int
ForceConstantMatrix::mem_layout( int ir, int mu2, int mu1) const
{
	assert( (ir >= 0) && (ir < this->get_num_R()) );
	int consqIndex = (ir*numModes_+mu2)*numModes_+mu1;
	assert( consqIndex < int(data_.size()) );
	return consqIndex;
}

int
ForceConstantMatrix::get_num_modes() const
{
	return numModes_;
}

int
ForceConstantMatrix::get_num_R() const
{
	return RVectorDim_[0]*RVectorDim_[1]*RVectorDim_[2];
}

void
ForceConstantMatrix::fourier_transform_q(std::vector<double> const & qVect,
		std::vector<std::complex<double>> & data) const
{
	assert( qVect.size()%3 == 0 );
	int nq = int(qVect.size())/3;
	int nM = this->get_num_modes();
	int nAtoms = nM/3;

	data.assign( nq*nM*nM , 0.0 );

	for ( int iq = 0 ; iq < nq ; ++iq)
	{
		std::vector< std::complex<double> > atomPosPhase(nAtoms);
		for ( int ia = 0 ; ia <nAtoms; ++ia )
			atomPosPhase[ia] = std::exp( std::complex<double>(0,
					2*M_PI*(qVect[iq*3]*tau_[0]+qVect[iq*3+1]*tau_[1]+qVect[iq*3+2]*tau_[2])) );

		for ( int iRz = 0 ; iRz < RVectorDim_[2]; ++iRz )
			for ( int iRy = 0 ; iRy < RVectorDim_[1]; ++iRy )
				for ( int iRx = 0 ; iRx < RVectorDim_[0]; ++iRx )
				{
					int R[3] = { iRx < RVectorDim_[0]/2 ? iRx : iRx - RVectorDim_[0],
									iRy < RVectorDim_[1]/2 ? iRy : iRy - RVectorDim_[1],
									iRz < RVectorDim_[2]/2 ? iRz : iRz - RVectorDim_[2] };
					double multi = RSupercellMultiplicity_[this->RVectorLayout(iRz,iRy,iRx)];
					std::complex<double> phase = std::exp( std::complex<double>(0,
							2*M_PI*(qVect[iq*3]*R[0]+qVect[iq*3+1]*R[1]+qVect[iq*3+2]*R[2])) ) / multi;

					//only compute upper half + diagonal
					for ( int mu1 = 0 ; mu1 < nM; ++mu1 )
						for ( int mu2 = mu1 ; mu2 < nM; ++mu2 )
							data[iq*nM*nM + mu1*nM+mu2] += atomPosPhase[mu1/3]*std::conj(atomPosPhase[mu2/3])*
														phase*data_[this->mem_layout(iRz,iRy,iRx,mu2,mu1)];
				}

		//reconstruct lower half as the conjugate
		for ( int mu1 = 0 ; mu1 < nM; ++mu1 )
			for ( int mu2 = 0 ; mu2 < mu1; ++mu2 )
				data[iq*nM*nM + mu1*nM+mu2] = std::conj(data[iq*nM*nM + mu2*nM+mu1]);
	}
}


void
ForceConstantMatrix::shift_atoms_into_map( std::vector<LatticeStructure::Atom> const & atoms,
		LatticeStructure::Atom const & site,
		std::map<LatticeStructure::Atom,int> & shiftedLattice) const
{
	shiftedLattice.clear();
	auto nCenter = site.get_position();
	auto hint = shiftedLattice.end();
	for ( int ia = 0 ; ia < atoms.size(); ++ia )
	{
		auto a = atoms[ia];
		auto p = a.get_position();
		for ( int i = 0 ; i < 3 ; ++i)
			p[i] -= nCenter[i];
		a.set_position( p );
		//optimal insertion if atoms is already sorted according to the internal atom comparison
		hint = shiftedLattice.insert(hint, std::move( std::make_pair( std::move(a), ia ) ) );
	}
}

void
ForceConstantMatrix::transform_map(
		std::map<LatticeStructure::Atom,int> & shiftedLattice,
		LatticeStructure::Symmetry::SymmetryOperation const & symOp,
		std::vector<int> & rotAtomsMap) const
{
	rotAtomsMap.resize( shiftedLattice.size() );
	for ( auto pair : shiftedLattice )
	{
		auto a = pair.first;
		a.transform(symOp);
		auto it = shiftedLattice.find( a );
		if ( it == shiftedLattice.end() )
			throw std::logic_error("The set of atoms is not closed under symmetry operations which can't be.");
		rotAtomsMap[pair.second] = it->second;
	}
}

void
ForceConstantMatrix::set_supercell_dim(
		LatticeStructure::UnitCell const & unitCell,
		LatticeStructure::UnitCell const & superCell)
{
	//Locate the unit cell in the supercell
	auto A = unitCell.get_lattice().get_latticeMatrix();
	for ( auto &aij : A )
		aij *= unitCell.get_alat();

	auto As = superCell.get_lattice().get_latticeMatrix();
	for ( auto &aij : As )
		aij *= superCell.get_alat();

	auto slice = [] (std::vector<double> const & A, int i) {
			return std::vector<double>({A[0*3+i],A[1*3+i],A[2*3+i]});
		};
	auto dot_p = [] (std::vector<double> const & a, std::vector<double> const & b) {
		return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
	};

	double scaleX  = dot_p(slice(A,0),slice(As,0))/dot_p(slice(A,0),slice(A,0));
	double scaleY  = dot_p(slice(A,1),slice(As,1))/dot_p(slice(A,1),slice(A,1));
	double scaleZ  = dot_p(slice(A,2),slice(As,2))/dot_p(slice(A,2),slice(A,2));

	supercellDim_ = std::vector<int> {
								int(std::floor( scaleX + 0.5 )),
								int(std::floor( scaleY + 0.5 )),
								int(std::floor( scaleZ + 0.5 )) };
	assert( (std::abs(scaleX - supercellDim_[0]) < superCell.get_symmetry().get_symmetry_prec()) &&
			(std::abs(scaleY - supercellDim_[1]) < superCell.get_symmetry().get_symmetry_prec()) &&
			(std::abs(scaleZ - supercellDim_[2]) < superCell.get_symmetry().get_symmetry_prec()) );

	RVectorDim_ = supercellDim_;
	for ( auto &Ri : RVectorDim_ )
		++Ri;
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

} /* namespace PhononStructure */
} /* namespace elephon */
