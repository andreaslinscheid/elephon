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
#include <map>
#include <cmath>
#include <stdexcept>
#include <string>
#include <set>

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
	symmetry_ = std::move(sym);
	//synchronize symmetries and lattice
	this->set_symmetry_to_lattice( symmetry_);
	this->generate_site_symmetries(atoms_,symmetry_,siteSymmetries_);
}

void
UnitCell::set_symmetry_to_lattice(LatticeStructure::Symmetry & symmetry) const
{
	symmetry.reset_lattice( lattice_ );
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
	for ( int iscx = 0 ; iscx < scx ; ++iscx)
		for ( int iscy = 0 ; iscy < scy ; ++iscy)
			for ( int iscz = 0 ; iscz < scz ; ++iscz)
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

	//The best strategy to save numerical effort is to reduce the number of irreducible displacements.
	//Thus, one should take displacements that are mapped to linear independent displacements by a site symmetry.
	//On the other hand, in order to make the calculation of the given displacement more efficient, the
	//overall symmetry group should not be affected.
	std::vector<AtomDisplacement> nonRotated;
	nonRotated.reserve( (symmetricDisplacements ? 6 : 3 ) );

	//Insert displacement along x (and -x)
	V vx = {1,0,0};
	this->add_displacement( vx, pos, symmetricDisplacements, name, delta, displMagn, nonRotated);

	//Insert displacement along y (and -y)
	V vy = {0,1,0};
	this->add_displacement( vy, pos, symmetricDisplacements, name, delta, displMagn, nonRotated);

	//Insert displacement along z (and -z)
	V vz = {0,0,1};
	this->add_displacement( vz, pos, symmetricDisplacements, name, delta, displMagn, nonRotated);

	//expand to all reducible displacements at this site
	std::vector<AtomDisplacement> rotated;
	rotated.reserve( nonRotated.size() * siteSymmetry.get_num_symmetries() );
	for ( int issym = 0 ; issym < siteSymmetry.get_num_symmetries(); ++issym )
		for ( auto ad : nonRotated)
		{
			ad.transform_direction( siteSymmetry.get_sym_op(issym) );
			rotated.push_back( std::move(ad) );
		}

	//remove equivalent displacements and find the irreducible set
	std::set<AtomDisplacement> reducibleSet( rotated.begin(), rotated.end() );
	reducible.assign( reducibleSet.begin(), reducibleSet.end() );
	SymmetryReduction<AtomDisplacement>( siteSymmetry,
			reducible,
			irreducible,
			redToIrredDispl, symRedToIrredDispl,
			irredToRedDispl, symIrredToRedDispl);
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
	lattice_.direct_to_cartesian(direction);
	double norm = magn(direction);
	for ( int i = 0; i < 3; ++i )
		direction[i] /= norm;

	//Create an atom displacement object
	addtothis.push_back( AtomDisplacement(atomName, magnInAngstroem, position, direction, gridPrec, (not symmetricDispl) ) );

	if (symmetricDispl)
	{
		//add the inverse direction, too
		for ( int i = 0; i < 3; ++i )
			direction[i] *= -1;
		addtothis.push_back( AtomDisplacement(atomName, magnInAngstroem, position, direction, gridPrec, (not symmetricDispl) ) );
	}
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

} /* namespace LatticeStructure */
} /* namespace elephon */
