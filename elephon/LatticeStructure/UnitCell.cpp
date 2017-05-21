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
	this->reduce_symmetry_to_lattice();
}

void UnitCell::reduce_symmetry_to_lattice()
{
	const double d = symmetry_.get_symmetry_prec();

	typedef std::vector<double> V;

	//Define an ordering for 3D points
	auto cmp = [d] (double a, double b) {
		if ( std::fabs( a-b ) < d )
			return 0;
		return ( a < b ? 1 : -1 );
	};

	auto compare_location = [cmp] (V const& l1, V const& l2) {
		for (int xi = 0 ; xi < 3; ++xi)
			if ( cmp(l1[xi],l2[xi]) )
				return l1[xi] < l2[xi];
		return false;
	};

	std::map< V,std::string, decltype(compare_location) > atomsSet(compare_location);
	for ( auto a : atoms_ )
	{
		auto ret = atomsSet.insert( std::make_pair(a.get_position(),a.get_kind()) );
		if ( not ret.second )
			throw std::logic_error( "Atom positions are not distinct as to symmetry precision" );
	}

	//Now we check if the atomsSet it mapped to itself for each symmetry operation
	std::vector<int> dropSym;
	for (int isym = 0 ; isym < symmetry_.get_num_symmetries(); ++isym)
		for ( auto a : atoms_ )
		{
			std::vector<double> position = a.get_position();
			symmetry_.apply(isym,position);
			auto ret = atomsSet.find(position);
			if ( ret == atomsSet.end() )
			{
				dropSym.push_back(isym);
				break;
			}
			if ( a.get_kind().compare(ret->second) != 0 )
			{
				dropSym.push_back(isym);
				break;
			}
		}

	//Check if we need to update the symmetry set
	if ( not dropSym.empty() )
		symmetry_.symmetry_reduction( dropSym );
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

void UnitCell::generate_displacements( double displMagn,
		bool symmetricDisplacements,
		std::vector<AtomDisplacement> & reducibleDisplacements,
		std::vector<AtomDisplacement> & irreducibleDisplacements,
		std::vector<int> & redToIrred,
		std::vector<int> & symRedToIrred,
		std::vector< std::vector<int> > & irredToRed,
		std::vector< std::vector<int> > & symIrredToRed ) const
{
	typedef std::vector<double> V;
	const double gridPrc = symmetry_.get_symmetry_prec();
	const double a0 = displMagn/lattice_.get_alat();
	const bool treatDirSymmetric = not symmetricDisplacements;
	auto magn = [] (V const& v) { return std::sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] ); };

	//First we generate all possible displacements (reducible)
	std::set<AtomDisplacement> reducibleSet;
	for ( auto a : atoms_)
	{
		//Insert displacement along x (and -x)
		if ( not a.get_movement_fixed()[0] )
		{
			V vx = {1,0,0};
			lattice_.direct_to_cartesian(vx);
			double dx = a0/magn(vx);
			reducibleSet.insert( AtomDisplacement(a.get_kind(), dx, a.get_position(),V({1,0,0}), gridPrc, treatDirSymmetric ) );
			if (symmetricDisplacements)
				reducibleSet.insert( AtomDisplacement(a.get_kind(), dx, a.get_position(),V({-1,0,0}), gridPrc, treatDirSymmetric ) );
		}

		//Insert displacement along y (and -y)
		if ( not a.get_movement_fixed()[1] )
		{
			V vy = {0,1,0};
			lattice_.direct_to_cartesian(vy);
			double dy = a0/magn(vy);
			reducibleSet.insert( AtomDisplacement(a.get_kind(), dy, a.get_position(),V({0,1,0}), gridPrc, treatDirSymmetric ) );
			if (symmetricDisplacements)
				reducibleSet.insert( AtomDisplacement(a.get_kind(), dy, a.get_position(),V({0,-1,0}), gridPrc, treatDirSymmetric ) );
		}

		//Insert displacement along z (and -z)
		if ( not a.get_movement_fixed()[2] )
		{
			V vz = {0,0,1};
			lattice_.direct_to_cartesian(vz);
			double dz = a0/magn(vz);
			reducibleSet.insert( AtomDisplacement(a.get_kind(), dz, a.get_position(),V({0,0,1}), gridPrc, treatDirSymmetric ) );
			if (symmetricDisplacements)
				reducibleSet.insert( AtomDisplacement(a.get_kind(), dz, a.get_position(),V({0,0,-1}), gridPrc, treatDirSymmetric ) );
		}
	}

	reducibleDisplacements = std::vector<AtomDisplacement>(
			reducibleSet.begin(),reducibleSet.end());

	SymmetryReduction<AtomDisplacement>(
			symmetry_,
			reducibleDisplacements,  irreducibleDisplacements,
			redToIrred, symRedToIrred,
			irredToRed, symIrredToRed);
}

void UnitCell::displace_atom( AtomDisplacement const& displ )
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

		//located atom
		std::vector<double> newPos = a.get_position();
		for (int i = 0 ; i < 3; ++i )
			newPos[i] += displ.generate_movement()[i];
		a.set_position( newPos );

		this->reduce_symmetry_to_lattice();
		return;
	}
	throw std::logic_error("Displacement does not match any atom");
}

std::vector<double> UnitCell::get_lattice_matrix() const
{
	return lattice_.get_latticeMatrix();
}

double UnitCell::get_alat() const
{
	return lattice_.get_alat();
}

} /* namespace LatticeStructure */
} /* namespace elephon */
