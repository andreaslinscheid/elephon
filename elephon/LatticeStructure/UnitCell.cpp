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
#include "symmetry/atom_transform_map.h"
#include "LatticeStructure/AtomSymmetryConnection.h"
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

	atomSymmetryModule_->initialize(atoms_, symmetry_);
}

void
UnitCell::set_symmetry_to_lattice(LatticeStructure::Symmetry & symmetry) const
{
	//This is not a trivial process because the symmetries are given in the lattice basis
	// and thus depend on the lattice. We first remove incompatible symmetries in the
	// lattice representation and then reset the cartesian representation.
	// TODO clean up and disentangle the symmetries. In principle the lattice could change the symmetries,
	//		too. We do not allow nor check for this at this point.
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
	std::vector<int> scale{scx,scy,scz};

	//Here we build the new lattice module
	auto newLattice = lattice_.build_supercell(scale);

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
UnitCell::compute_supercell_dim(
		std::shared_ptr<const UnitCell> superCell,
		std::vector<int> & supercellDim ) const
{
	//Locate the unit cell in the supercell
	auto dot_p = [] (std::vector<double> const & a, std::vector<double> const & b) {
		return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
	};

	auto a1 = lattice_.get_lattice_vector(0);
	auto a2 = lattice_.get_lattice_vector(1);
	auto a3 = lattice_.get_lattice_vector(2);

	double scaleX  = superCell->get_lattice().get_alat()/lattice_.get_alat()
			*dot_p(a1,superCell->get_lattice().get_lattice_vector(0))/dot_p(a1,a1);
	double scaleY  = superCell->get_lattice().get_alat()/lattice_.get_alat()
			*dot_p(a2,superCell->get_lattice().get_lattice_vector(1))/dot_p(a2,a2);
	double scaleZ  = superCell->get_lattice().get_alat()/lattice_.get_alat()
			*dot_p(a3,superCell->get_lattice().get_lattice_vector(2))/dot_p(a3,a3);

	supercellDim = std::vector<int> {
								int(std::floor( scaleX + 0.5 )),
								int(std::floor( scaleY + 0.5 )),
								int(std::floor( scaleZ + 0.5 )) };
	assert( (std::abs(scaleX - supercellDim[0]) < superCell->get_symmetry().get_symmetry_prec()) &&
			(std::abs(scaleY - supercellDim[1]) < superCell->get_symmetry().get_symmetry_prec()) &&
			(std::abs(scaleZ - supercellDim[2]) < superCell->get_symmetry().get_symmetry_prec()) );
}

std::shared_ptr<const LatticeStructure::AtomSymmetryConnection>
UnitCell::get_atom_symmetry() const
{
	return atomSymmetryModule_;
}

} /* namespace LatticeStructure */
} /* namespace elephon */
