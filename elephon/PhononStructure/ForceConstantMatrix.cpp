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

	assert( forces.size() == irredDispl.size() );

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

		std::vector<LatticeStructure::AtomDisplacement> irredNotUsed,reducible;
		bool symmetricDispl = not irredDispl[iForceField].is_plus_minus_displ_equivalent();
		std::vector<int> redToIrredDispl, symRedToIrredDispl;
		std::vector< std::vector<int> > irredToRedDispl, symIrredToRedDispl;
		unitCell.get_site_displacements( irredAtoms[ia],
				symmetricDispl, sSym,
				1, // not used here
				irredNotUsed,reducible,
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
		auto p = irredAtoms[ia].get_position();
		for ( int i = 0 ; i < 3 ; ++i)
			p[i] /= double(supercellDim_[i]);

		//For reducible displacements in the star of an irreducible one, this means, we need to
		//	1) identify position g_\alpha^-1 * i, where i is a supercell atom index and g_\alpha is a rotation
		//		around the current site index.
		//	2) rotate the 3-vector with the forces by g_\alpha
		//For all reducible displacement at this site, build the transpose of the matrix F
		//NOTE that the inverse symmetry operation has to rotate first and displace later.
		std::vector< std::vector<int> > rotMapsSiteSymmetry;
		this->transform_map( p, superCell.get_atoms_list(), sSym, rotMapsSiteSymmetry);
//for ( int isym = 0 ; isym <sSym.get_num_symmetries(); ++isym)
//{
//	for ( int iaS = 0 ; iaS < naSC ; ++iaS )
//		std::cout << ' ' << rotMapsSiteSymmetry[isym][iaS];
//	std::cout << std::endl;
//}

		std::vector<double> matrixForceConstants(3*3);
		std::vector<double> F( 3*nRedDispl );
		for ( int iaS = 0 ; iaS < naSC ; ++iaS )
		{
			for ( int idir = 0 ; idir < irredToRedDispl.size(); ++idir)
			{
				for ( int istar = 0 ; istar < symIrredToRedDispl[idir].size(); ++istar)
				{
					int id = irredToRedDispl[idir][istar];
					int iForceAtom = rotMapsSiteSymmetry[symIrredToRedDispl[idir][istar]][iaS];
					std::copy( &forces[iForceField+idir][iForceAtom*3],
							&forces[iForceField+idir][iForceAtom*3]+3,
							&F[id*3] );
					sSym.rotate_cartesian( symIrredToRedDispl[idir][istar],
							F.begin()+id*3, F.begin()+(id+1)*3 );
				}
			}

//std::cout <<'\n'<< iaS<< std::endl;
//for ( int i = 0 ; i < U.size()/3 ; ++i)
//{
//	for ( int j = 0 ; j < 3 ; ++j)
//		std::cout << U[i*3+j] << '\t';
//	std::cout << std::endl;
//}
//std::cout << std::endl;
//for ( int i = 0 ; i < F.size()/3 ; ++i)
//{
//	for ( int j = 0 ; j < 3 ; ++j)
//		std::cout << F[i*3+j] << '\t';
//	std::cout << std::endl;
//}
			//multiplying the Moore-Penrose pseudo inverse, of the displacements,
			//we obtain the least square fit of transpose of the force constant matrix
			linAlg.matrix_matrix_prod( piU, F, matrixForceConstants, 3, 3 );

			//transpose back
			auto localMat = irreducibleMatrixOfForceConstants.data()+(ia*naSC+iaS)*9;
//std::cout << ia << '\t' << iaS<< std::endl;
			for ( int xi = 0 ; xi < 3; ++xi)
			{
				for ( int xj = 0 ; xj < 3; ++xj)
				{
					localMat[xi*3+xj] = -matrixForceConstants[xj*3+xi];
//std::cout << '\t' << localMat[xi*3+xj];
				}
//std::cout << std::endl;
			}
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
	std::vector<double> matForceSlice( 9*naSC );
	for ( int irA = 0; irA < irredToRedAtoms.size(); ++irA )
	{
		//loop the star of the irreducible atom
		for ( int istar = 0; istar < symIrredToRedAtoms[irA].size(); ++istar )
		{
			//transform the preliminary matrix of force constants according
			//to the present symmetry operation which takes us from the irreducible
			//to the reducible atom
			int isym = symIrredToRedAtoms[irA][istar];

			//Re-shuffle the matrices according to the symmetry operation beta = S( alpha )
			//Note: at this point we need that the symmetry of the unit cell is the same as
			//		the one of the supercell
			for ( int iaSC = 0 ; iaSC < naSC ; ++iaSC )
			{
				auto a = superCell.get_atoms_list()[iaSC];
				a.transform(unitCell.get_symmetry().get_sym_op(isym));
				auto it = superCellLookup.find(a);
				if ( it == superCellLookup.end() )
					throw std::logic_error("Unable to locate rotates atom!");
				int iredA1 = it->second;
				std::copy(&irreducibleMatrixOfForceConstants[(irA*naSC+iaSC)*9],
						  &irreducibleMatrixOfForceConstants[(irA*naSC+iaSC)*9]+9,
						&matForceSlice[iredA1*9] );
			}

			//transform the set of 3x3 matrices for the set of atoms
			unitCell.get_symmetry().rotate_matrix_cartesian( isym,
					matForceSlice.begin(),
					matForceSlice.end());

			//Convert into the format of the dynamical matrix in real space
			for ( int iaSC = 0 ; iaSC < naSC ; ++iaSC )
			{
				int ia1 = atomSuperCellToPrimitve[iaSC].first;
				int iaDispl = irredToRedAtoms[irA][istar];
				std::vector<int> R = atomSuperCellToPrimitve[iaSC].second;
				for ( int i = 0 ; i < 3 ; ++i)
					R[i] = R[i] < 0 ? R[i] + supercellDim_[i] : R[i];
				int ir = this->RVectorLayout(R[2],R[1],R[0]);
//std::cout << ia1 << '\t' << iaDispl<< std::endl;
				for ( int xj = 0 ; xj < 3; ++xj)
				{
					for ( int xi = 0 ; xi < 3; ++xi)
					{
						data_[this->mem_layout(ir,ia1*3+xj,iaDispl*3+xi)] = matForceSlice[iaSC*9+xj*3+xi];
//std::cout << '\t' <<matForceSlice[iaSC*9+xj*3+xi];
					}
//std::cout << std::endl;
				}
			}
		}
	}

	int nAUC = unitCell.get_atoms_list().size();
	tau_.clear();
	tau_.resize( nAUC*naSC );
	for ( int ia = 0; ia < nAUC; ++ia )
	{
		auto atomsSC = superCell.get_atoms_list();
		auto shift = unitCell.get_atoms_list()[ia].get_position();
		for ( int i = 0 ; i < 3 ; ++i )
			shift[i] /= double(supercellDim_[i]);
		for ( auto &a : atomsSC)
		{
			auto pos = a.get_position();
			for ( int i = 0 ; i < 3 ; ++i )
				pos[i] -= shift[i];
			a.set_position(pos);
		}
		//atomsSC is now a supercell centered around atom ia

		//iff an atoms is close to the boarder we have to place a copy
		//on the respective other one to preserve the local symmetry [see PRL]
		for ( int ib = 0 ; ib < naSC; ++ib)
		{
			auto pos = atomsSC[ib].get_position();
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
			tau_[ia*naSC+ib].resize(vects.size()*3);
			for ( int iR = 0 ; iR < vects.size(); ++iR )
				for ( int xi = 0 ; xi < 3; ++xi)
					tau_[ia*naSC+ib][iR*3+xi] = vects[iR][xi]*supercellDim_[xi];
		}
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
ForceConstantMatrix::operator() (int Rz, int Ry, int Rx, int mu2, int mu1) const
{
	return data_[ this->mem_layout(Rz,Ry,Rx,mu2,mu1) ];
}

int
ForceConstantMatrix::RVectorLayout(int iRz, int iRy, int iRx ) const
{
	assert( (iRx >= 0) && (iRx < supercellDim_[0]) );
	assert( (iRy >= 0) && (iRy < supercellDim_[1]) );
	assert( (iRz >= 0) && (iRz < supercellDim_[2]) );
	return (iRz*supercellDim_[1]+iRy)*supercellDim_[0]+iRx;
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
	return supercellDim_[0]*supercellDim_[1]*supercellDim_[2];
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

	int nAUC = nM/3;
	int nR = this->get_num_R();

	std::vector< std::complex<double> > phases(nR*nAUC*nAUC);
	for ( int iq = 0 ; iq < nq ; ++iq)
	{
		for ( int ia = 0 ; ia < nAUC; ++ia )
			for ( int ib = 0 ; ib < nAUC; ++ib )
				for ( int ir = 0 ; ir < nR; ++ir )
				{
					int multi = tau_[(ia*nAUC+ib)*nR+ir].size()/3;
					phases[(ia*nAUC+ib)*nR+ir] = 0;
					for ( int im = 0 ; im < multi; ++im )
						phases[(ia*nAUC+ib)*nR+ir] += std::exp( std::complex<double>(0,
							-2*M_PI*(qVect[iq*3+0]*tau_[(ia*nAUC+ib)*nR+ir][im*3+0]
							        +qVect[iq*3+1]*tau_[(ia*nAUC+ib)*nR+ir][im*3+1]
									+qVect[iq*3+2]*tau_[(ia*nAUC+ib)*nR+ir][im*3+2])) ) / double(multi);
				}

		//only compute upper half + diagonal
		for ( int mu1 = 0 ; mu1 < nM; ++mu1 )
			for ( int mu2 = 0 ; mu2 < nM; ++mu2 )
			{
				int ia = mu1/3;
				int ib = mu2/3;

				for ( int ir = 0 ; ir < nR; ++ir )
				{

					data[iq*nM*nM + mu1*nM+mu2] += phases[(ia*nAUC+ib)*nR+ir]*data_[this->mem_layout(ir,mu1,mu2)];

				}
			}

//		//reconstruct lower half as the conjugate
//		for ( int mu1 = 0 ; mu1 < nM; ++mu1 )
//			for ( int mu2 = mu1+1 ; mu2 < nM; ++mu2 )
//				data[iq*nM*nM + mu2*nM+mu1] = std::conj(data[iq*nM*nM + mu1*nM+mu2]);
	}
}

void
ForceConstantMatrix::transform_map(
		std::vector<double> const & shift,
		std::vector<LatticeStructure::Atom> atoms,
		LatticeStructure::Symmetry const & siteSymmetry,
		std::vector< std::vector<int> > & rotAtomsMap) const
{
	int iA = atoms.size();
	int iS = siteSymmetry.get_num_symmetries();
	//We shift and create a lookup for atoms, then we apply the rotation and
	//try to discover the transformed position in the set
	std::map<LatticeStructure::Atom,int> loopup;
	for ( int i = 0 ; i < iA ; ++i)
	{
		auto pr = atoms[i].get_position();
		for ( int xi = 0 ; xi < 3; ++xi)
			pr[xi] -= shift[xi];
		atoms[i].set_position(pr);
		loopup.insert( std::move(std::make_pair(atoms[i],i)) );
	}

	rotAtomsMap = std::vector< std::vector<int> >(iS, std::vector<int>(iA) );
	for ( int isym = 0 ; isym < iS; ++isym)
	{
		//rotate all atoms
		auto rotAtoms = atoms;
		for ( auto &a : rotAtoms )
			a.transform(siteSymmetry.get_sym_op(isym));

		for ( int i = 0 ; i < iA ; ++i )
		{
			auto it = loopup.find( rotAtoms[i] );
			if ( it == loopup.end() )
				throw std::logic_error("The set of atoms is not closed "
						"under symmetry operations which can't be.");
			//rotAtomsMap tells for a given atom where it ends up after application of the inverse symmetry
			rotAtomsMap[isym][it->second] = i;
		}
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
