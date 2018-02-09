/*	This file Atom.cpp is part of elephon.
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

#include "LatticeStructure/Atom.h"
#include "LatticeStructure/AtomDisplacement.h"
#include <cmath>
#include <assert.h>
#include <stdexcept>

namespace elephon
{
namespace LatticeStructure
{

Atom::Atom() { };

Atom::Atom(
		double mass,
		std::string kind,
		std::vector<double> pos,
		std::vector<bool> frozen,
		double gridPrec)
		: mass_(mass), kind_(std::move(kind)),pos_(std::move(pos)),frozen_(std::move(frozen)), prec_(gridPrec)
{
	assert(mass_>0);
	assert(not kind_.empty());
	assert(pos_.size() == 3);
	assert(frozen_.size() == 3);
	this->map_pos_back_1BZ();
};

std::string Atom::get_kind() const
{
	return kind_;
};

double
Atom::get_mass() const
{
	if ( mass_ == 0 )
		throw std::logic_error("Trying to get atomic mass that was not set");
	return mass_;
}

double
Atom::get_position_precision() const
{
	return prec_;
}

std::vector<double> const & Atom::get_position() const
{
	assert(pos_.size() == 3);
	return pos_;
};

void Atom::set_position(std::vector<double> newPos)
{
	pos_ = std::move(newPos);
	this->map_pos_back_1BZ();
}

std::vector<bool> Atom::get_movement_fixed() const
{
	assert(frozen_.size() == 3);
	return frozen_;
};

void
Atom::transform( symmetry::SymmetryOperation const & sop  )
{
	sop.apply(pos_,true);
}

void
Atom::apply_displacement(AtomDisplacement const & displ)
{
	assert(pos_.size()==3);
	if(not ((std::abs(pos_[0]-displ.get_position()[0])<prec_)&&
			(std::abs(pos_[1]-displ.get_position()[1])<prec_)&&
			(std::abs(pos_[2]-displ.get_position()[2])<prec_) ) )
		throw std::logic_error("Displacement applied to atom not referenced by the displacement");

	for (int i = 0 ; i < 3 ; ++i)
		pos_[i] += displ.get_direction()[i]*displ.get_magnitude();
	this->map_pos_back_1BZ();
}

void
Atom::map_pos_back_1BZ()
{
	assert(pos_.size()==3);
	for ( auto &xi : pos_ )
		xi -= std::floor(xi + 0.5);
}

bool operator< (Atom const & a1, Atom const & a2)
{
	assert(a1.prec_ == a2.prec_);
	for ( int i = 3; i-- ; )
		if ( std::abs(a1.pos_[i] - a2.pos_[i]) >= a1.prec_ )
			return a1.pos_[i] < a2.pos_[i];
	return false;
}

bool operator== (Atom const & a1, Atom const & a2)
{
	bool equalPosition = (!(a1<a2))&&(!(a2<a1));
	// position is equal but kind is differen => different atoms (? throw here - this makes probably no sense ?)
	if( equalPosition and (a1.get_kind().compare(a2.get_kind()) != 0 ))
		return false;
	return equalPosition;
};


} /* namespace LatticeStructure */
} /* namespace elephon */
