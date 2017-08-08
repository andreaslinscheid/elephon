/*	This file UnitCell.cpp is part of elephon.
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

#include "UnitCell.h"
#include "SymmetryReduction.h"
#include "Algorithms/LinearAlgebraInterface.h"
#include <map>
#include <cmath>
#include <stdexcept>
#include <string>
#include <set>
#include <iostream>
#include <algorithm>

namespace elephon
{
namespace LatticeStructure
{

void UnitCell::initialize(
		std::vector<LatticeStructure::Atom> atoms,
		LatticeStructure::LatticeModule lattice,
		LatticeStructure::Symmetry sym)
{
	atoms_ = std::move(atoms);
	lattice_ = std::move(lattice);
	assert( not sym.is_reci() );
	symmetry_ = std::move(sym);
	//synchronize symmetries and lattice
	this->set_symmetry_to_lattice( symmetry_);
	this->generate_site_symmetries(atoms_,symmetry_,siteSymmetries_);
}

void
UnitCell::set_symmetry_to_lattice(LatticeStructure::Symmetry & symmetry) const
{
	//This is not a trivial process because the symmetries are given in the lattice basis
	// and thus depend on the lattice. We first remove incompatible symmetries in the
	// lattice representation and then reset the cartesian representation.
	// TODO clean up and disentangle the symmetries. In principle the lattice could change the symmetries,
	//		too. We do not allow nor check for this at this point.
	const double d = symmetry.get_symmetry_prec();

	typedef std::vector<double> V;

	std::set< Atom > atomsSet;
	for ( auto a : atoms_ )
	{
		auto ret = atomsSet.insert( std::move(a) );
		if ( not ret.second )
			throw std::logic_error( "Atom positions are not distinct as to symmetry precision" );
	}

	//Now we check if the atomsSet it mapped to itself for each symmetry operation
	std::vector<int> dropSym;
	for (int isym = 0 ; isym < symmetry.get_num_symmetries(); ++isym)
		for ( auto a : atoms_ )
		{
			a.transform( symmetry.get_sym_op( isym ) );
			auto ret = atomsSet.find(a);
			if ( ret == atomsSet.end() )
			{
				dropSym.push_back(isym);
				break;
			}
			if ( a.get_kind().compare(ret->get_kind()) != 0 )
			{
				dropSym.push_back(isym);
				break;
			}
		}

	//Check if we need to update the symmetry set
	if ( not dropSym.empty() )
		symmetry.symmetry_reduction( dropSym );
	symmetry.reset_lattice( lattice_ );
}

UnitCell UnitCell::build_supercell(int scx, int scy, int scz) const
{
	double scale[3] = {double(scx),double(scy),double(scz)};

	//Here we build the new lattice module
	std::vector<double> newLatticeMatrix = lattice_.get_latticeMatrix();
	for ( int i = 0 ; i < 3; ++i)
		for ( int j = 0 ; j < 3; ++j)
			newLatticeMatrix[i*3+j] *= scale[i]*lattice_.get_alat();
	LatticeStructure::LatticeModule newLattice;
	newLattice.initialize(newLatticeMatrix);

	//Here we build the new set of Atoms
	std::vector<LatticeStructure::Atom> newAtoms;
	newAtoms.reserve( atoms_.size()*scx*scy*scz );
	for ( int iscz = 0 ; iscz < scz ; ++iscz)
		for ( int iscy = 0 ; iscy < scy ; ++iscy)
			for ( int iscx = 0 ; iscx < scx ; ++iscx)
			{
				std::vector<double> baseVect =
					{ double(iscx)/double(scx),double(iscy)/double(scy),double(iscz)/double(scz) };

				for ( auto a : atoms_)
				{
					auto newPos = a.get_position( );
					for ( int i = 0 ; i < 3; ++i )
						newPos[i] = newPos[i]/scale[i] + baseVect[i];
					a.set_position(newPos);
					newAtoms.push_back( std::move(a) );
				}
			}

	//Possible symmetry reduction will take place automatically
	UnitCell supercell;
	supercell.initialize(newAtoms,newLattice,symmetry_);
	return supercell;
}

std::vector<LatticeStructure::Atom> const &
UnitCell::get_atoms_list() const
{
	return atoms_;
}

void
UnitCell::generate_displacements( double displMagn,
		bool symmetricDisplacements,
		std::vector<AtomDisplacement> & irreducibleDisplacements) const
{
	irreducibleDisplacements.clear();

	//find equivalent atoms
	std::vector<Atom> irredAtoms;
	std::vector<int> redToIrredAtoms;
	std::vector<int> symRedToIrredAtoms;
	std::vector< std::vector<int> > irredToRedAtoms;
	std::vector< std::vector<int> > symIrredToRedAtoms;
	SymmetryReduction<Atom>(
			symmetry_,
			this->get_atoms_list(),  irredAtoms,
			redToIrredAtoms, symRedToIrredAtoms,
			irredToRedAtoms, symIrredToRedAtoms);

	std::set<AtomDisplacement> reducibleSet;
	for ( int ia = 0 ; ia < irredAtoms.size(); ++ia)
	{
		std::vector<AtomDisplacement> irreducibleThisAtom, reducibleThisAtom;
		std::vector<int> redToIrred, symRedToIrred;
		std::vector< std::vector<int> > irredToRed, symIrredToRed;

		this->get_site_displacements( irredAtoms[ia], symmetricDisplacements, siteSymmetries_[ia], displMagn,
				irreducibleThisAtom, reducibleThisAtom,
				 redToIrred, symRedToIrred,
				 irredToRed, symIrredToRed);

		irreducibleDisplacements.insert( irreducibleDisplacements.end(),
				irreducibleThisAtom.begin(), irreducibleThisAtom.end() );
	}
}

void
UnitCell::get_site_displacements(Atom const & atomicSite,
		bool symmetricDisplacements,
		LatticeStructure::Symmetry const & siteSymmetry,
		double displMagn,
		std::vector<AtomDisplacement> & irreducible,
		std::vector<AtomDisplacement> & reducible,
		std::vector<int> & redToIrredDispl,
		std::vector<int> & symRedToIrredDispl,
		std::vector< std::vector<int> > & irredToRedDispl,
		std::vector< std::vector<int> > & symIrredToRedDispl) const
{
	typedef std::vector<double> V;
	auto pos = atomicSite.get_position();
	auto name = atomicSite.get_kind();

	double delta = siteSymmetry.get_symmetry_prec();
	Algorithms::LinearAlgebraInterface linAlg;

	//The best strategy to save numerical effort is to reduce the number of irreducible displacements.
	//Thus, one should take displacements that are mapped to linear independent displacements by a site symmetry.
	//On the other hand, in order to make the calculation of the given displacement more efficient, the
	//overall symmetry group should not be affected.
	std::vector<AtomDisplacement> nonRotated;
	std::vector<AtomDisplacement> rotated;

	auto symmetry_expand = [&] () {
		rotated.clear();
		rotated.reserve( nonRotated.size() * siteSymmetry.get_num_symmetries() );
		irredToRedDispl.resize(nonRotated.size(), std::vector<int>( siteSymmetry.get_num_symmetries() ));
		symIrredToRedDispl.resize(nonRotated.size(), std::vector<int>( siteSymmetry.get_num_symmetries() ));
		redToIrredDispl.resize(nonRotated.size()* siteSymmetry.get_num_symmetries());
		symRedToIrredDispl.resize(nonRotated.size()* siteSymmetry.get_num_symmetries());
		for ( int issym = 0 ; issym < siteSymmetry.get_num_symmetries(); ++issym )
		{
			for (int inr = 0 ; inr < nonRotated.size(); ++inr)
			{
				irredToRedDispl[inr][issym] = rotated.size();
				symIrredToRedDispl[inr][issym] = issym;
				redToIrredDispl[rotated.size()] = inr;
				symRedToIrredDispl[rotated.size()] = siteSymmetry.get_index_inverse(issym);
				auto ad = nonRotated[inr];
				ad.transform_direction( siteSymmetry.get_sym_op(issym) );
				rotated.push_back( std::move(ad) );
			}
		}
	};

	//Insert displacement along lattice x
	V vx = {1,0,0};
	lattice_.direct_to_cartesian_angstroem(vx);
	nonRotated.push_back( AtomDisplacement(name, displMagn, pos, vx, delta, (not symmetricDisplacements) ) );

	//expand to all reducible displacements at this site
	symmetry_expand();

	//Add the minus displacement if it is not already in the set
	V mvx = {-1,0,0};
	lattice_.direct_to_cartesian_angstroem(mvx);
	AtomDisplacement minusDx(name, displMagn, pos, mvx, delta, (not symmetricDisplacements) );
	auto it = std::find( rotated.begin(), rotated.end(), minusDx);
	if ( (it == rotated.end()) and (symmetricDisplacements) )
	{
		nonRotated.push_back( std::move(minusDx) );
		symmetry_expand();
	}

	V redDisplTmp;
	redDisplTmp.resize(3*rotated.size());
	for ( int i = 0 ; i < rotated.size(); ++i)
		for ( int j = 0 ; j < 3; ++j)
		redDisplTmp[i*3+j] = rotated[i].get_direction()[j];

	//Check if the resulting set of vectors spans full 3D space.
	//Very small overlaps (<0.01) are considered null
	int nullDim;
	V nullSpace;
	linAlg.null_space(redDisplTmp,rotated.size(),3,nullDim,nullSpace, 0.01);

	auto normalize = [] (V & v ){
		assert(v.size() == 3 );
		double norm = std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
		for ( auto &a : v)
			a /= norm;
	};

	auto check_overlap_is_small = [] ( V const & nullSpace, int nullDim, V const & trialVector) {
		const double small = 0.01;
		assert(nullSpace.size() == 3*nullDim );
		assert(trialVector.size() == 3 );
		bool overlapIsZero = true;
		for ( int i = 0 ; i < nullDim ; ++i )
		{
			double overlap = 0;
			for ( int j = 0 ; j < 3 ; ++j )
				overlap += nullSpace[i*3+j]*trialVector[j];
			overlapIsZero = overlapIsZero and (overlap < small);
		}
		return overlapIsZero;
	};

	if ( nullDim != 0 )
	{
		//x + rotations do not span 3D space, insert displacement
		// along y (and -y) if its overlap with the nullSpace is not zero
		V vy = {0,1,0};
		lattice_.direct_to_cartesian_angstroem(vy);
		if ( not check_overlap_is_small(nullSpace, nullDim, vy) )
			nonRotated.push_back( AtomDisplacement(name, displMagn, pos, vy, delta, (not symmetricDisplacements) ) );

		//nonRotated contains now x and y - expand by symmetry and check null space again
		symmetry_expand();

		V mvy = {0,-1,0};
		lattice_.direct_to_cartesian_angstroem(mvy);
		AtomDisplacement minusDy(name, displMagn, pos, mvy, delta, (not symmetricDisplacements) );
		auto it = std::find( rotated.begin(), rotated.end(), minusDy);
		if ( (it == rotated.end()) and (symmetricDisplacements) )
		{
			nonRotated.push_back( std::move(minusDy) );
			symmetry_expand();
		}

		redDisplTmp.resize(3*rotated.size());
		for ( int i = 0 ; i < rotated.size(); ++i)
			for ( int j = 0 ; j < 3; ++j)
			redDisplTmp[i*3+j] = rotated[i].get_direction()[j];
		linAlg.null_space(redDisplTmp,rotated.size(),3,nullDim,nullSpace, 0.01);
	}

	if ( nullDim != 0 )
	{
		//x, y + rotations do not span 3D space, insert displacement
		// along z (and -z) if its overlap with the nullSpace is not zero
		V vz = {0,0,1};
		lattice_.direct_to_cartesian_angstroem(vz);
		nonRotated.push_back( AtomDisplacement(name, displMagn, pos, vz, delta, (not symmetricDisplacements) ) );

		symmetry_expand();

		V mvz = {0,0,-1};
		lattice_.direct_to_cartesian_angstroem(mvz);
		AtomDisplacement minusDz(name, displMagn, pos, mvz, delta, (not symmetricDisplacements) );
		auto it = std::find( rotated.begin(), rotated.end(), minusDz);
		if ( (it == rotated.end()) and (symmetricDisplacements) )
		{
			nonRotated.push_back( std::move(minusDz) );
			symmetry_expand();
		}

		//In principle, if the algorithm work, the following is unnecessary but better check
		//nonRotated contains now x, y and z - expand by symmetry and check null space again
		redDisplTmp.resize(3*rotated.size());
		for ( int i = 0 ; i < rotated.size(); ++i)
			for ( int j = 0 ; j < 3; ++j)
			redDisplTmp[i*3+j] = rotated[i].get_direction()[j];
		linAlg.null_space(redDisplTmp,rotated.size(),3,nullDim,nullSpace, 0.01);
	}

	if ( nullDim != 0 )
		throw std::logic_error( "Displacement generation algorithm failed to span all 3D space!" );

	irreducible  = std::move(nonRotated);
	reducible = std::move(rotated);
}

void
UnitCell::add_displacement( std::vector<double> direction,
		std::vector<double> const & position,
		bool symmetricDispl,
		std::string const & atomName,
		double gridPrec,
		double magnInAngstroem,
		std::vector<AtomDisplacement> & addtothis) const
{
	auto magn = [] (std::vector<double> const& v) { return std::sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] ); };

	//create normalized cartesian vector for the displacement direction
	lattice_.direct_to_cartesian_angstroem(direction);
	double norm = magn(direction);
	for ( int i = 0; i < 3; ++i )
		direction[i] /= norm;

	//Create an atom displacement object
	addtothis.push_back( AtomDisplacement(atomName, magnInAngstroem, position, direction, gridPrec, (not symmetricDispl) ) );
}

void
UnitCell::displace_atom( AtomDisplacement const& displ )
{
	for ( auto &a : atoms_ )
	{
		if ( std::abs(a.get_position()[0]-displ.get_position()[0]) >= displ.get_prec() )
			continue;
		if ( std::abs(a.get_position()[1]-displ.get_position()[1]) >= displ.get_prec() )
			continue;
		if ( std::abs(a.get_position()[2]-displ.get_position()[2]) >= displ.get_prec() )
			continue;

		if ( a.get_kind().compare( displ.get_kind() ) != 0 )
			throw std::logic_error(" Displacement does not match the atom at this site ");

		//displace located atom
		std::vector<double> newPos = a.get_position();
		auto v = displ.get_direction();
		lattice_.cartesian_to_direct( v );
		for (int i = 0 ; i < 3; ++i )
			newPos[i] += v[i]*displ.get_magnitude()/lattice_.get_alat();
		a.set_position( newPos );

		this->set_symmetry_to_lattice(symmetry_);
		return;
	}
	throw std::logic_error("Displacement does not match any atom");
}

LatticeStructure::LatticeModule const &
UnitCell::get_lattice() const
{
	return lattice_;
}

double UnitCell::get_alat() const
{
	return lattice_.get_alat();
}

LatticeStructure::Symmetry const &
UnitCell::get_symmetry() const
{
	return symmetry_;
}

void
UnitCell::generate_site_symmetries(std::vector<LatticeStructure::Atom> const & atoms,
		LatticeStructure::Symmetry const & symmetry,
	std::vector<LatticeStructure::Symmetry> & siteSymmetries) const
{
	siteSymmetries.clear();
	siteSymmetries.reserve( atoms.size() );
	for ( int ia = 0 ; ia < atoms.size(); ++ia)
	{
		auto ss = symmetry;
		ss.small_group( atoms[ia].get_position() );
		siteSymmetries.push_back( std::move(ss) );
	}
}

LatticeStructure::Symmetry const &
UnitCell::get_site_symmetry(int atomIndex) const
{
	assert( (atomIndex >= 0) && (atomIndex < siteSymmetries_.size() ) );
	return siteSymmetries_[atomIndex];
}

void
UnitCell::generate_rotation_maps(std::vector<std::vector<int> > & rotationMap) const
{
	std::map< Atom, int > cmpMap;
	for ( int ia = 0 ; ia < atoms_.size() ; ++ia)
		cmpMap.insert( std::make_pair( atoms_[ia], ia ) );

	rotationMap.resize( symmetry_.get_num_symmetries() );
	for (int isym = 0 ; isym < symmetry_.get_num_symmetries(); ++isym)
	{
		rotationMap[isym].resize( atoms_.size() );
		for ( int ia = 0 ; ia < atoms_.size() ; ++ia)
		{
			auto a = atoms_[ia];
			a.transform( symmetry_.get_sym_op(isym) );
			auto ret = cmpMap.find(a);
			if ( ret == cmpMap.end() )
				throw std::logic_error("Cannot find rotated atom: set of atoms not closed under rotations!");
			rotationMap[isym][ia] = ret->second;
		}
	}
}

void
UnitCell::compute_supercell_dim(
		UnitCell const & superCell,
		std::vector<int> & supercellDim ) const
{
	//Locate the unit cell in the supercell
	auto dot_p = [] (std::vector<double> const & a, std::vector<double> const & b) {
		return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
	};

	auto a1 = lattice_.get_lattice_vector(0);
	auto a2 = lattice_.get_lattice_vector(1);
	auto a3 = lattice_.get_lattice_vector(2);

	double scaleX  = superCell.get_lattice().get_alat()/lattice_.get_alat()
			*dot_p(a1,superCell.get_lattice().get_lattice_vector(0))/dot_p(a1,a1);
	double scaleY  = superCell.get_lattice().get_alat()/lattice_.get_alat()
			*dot_p(a2,superCell.get_lattice().get_lattice_vector(1))/dot_p(a2,a2);
	double scaleZ  = superCell.get_lattice().get_alat()/lattice_.get_alat()
			*dot_p(a3,superCell.get_lattice().get_lattice_vector(2))/dot_p(a3,a3);

	supercellDim = std::vector<int> {
								int(std::floor( scaleX + 0.5 )),
								int(std::floor( scaleY + 0.5 )),
								int(std::floor( scaleZ + 0.5 )) };
	assert( (std::abs(scaleX - supercellDim[0]) < superCell.get_symmetry().get_symmetry_prec()) &&
			(std::abs(scaleY - supercellDim[1]) < superCell.get_symmetry().get_symmetry_prec()) &&
			(std::abs(scaleZ - supercellDim[2]) < superCell.get_symmetry().get_symmetry_prec()) );
}

} /* namespace LatticeStructure */
} /* namespace elephon */
