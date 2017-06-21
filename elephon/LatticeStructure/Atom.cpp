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

#include "Atom.h"
#include <cmath>
#include <assert.h>

namespace elephon
{
namespace LatticeStructure
{

Atom::Atom(std::string kind, std::vector<double> pos,
		std::vector<bool> frozen,
		double gridPrec)
		: kind_(std::move(kind)),pos_(std::move(pos)),frozen_(std::move(frozen)), prec_(gridPrec)
{
	for ( auto &xi : pos_ )
		xi -= std::floor(xi + 0.5);
};

std::string Atom::get_kind() const
{
	return kind_;
};

std::vector<double> const & Atom::get_position() const
{
	return pos_;
};

void Atom::set_position(std::vector<double> newPos)
{
	pos_ = std::move(newPos);
	for ( auto &xi : pos_ )
		xi -= std::floor(xi + 0.5);
}

std::vector<bool> Atom::get_movement_fixed() const
{
	return frozen_;
};

void
Atom::transform( Symmetry::SymmetryOperation const & sop  )
{
	sop.apply(pos_,true);
}

bool operator< (Atom const & a1, Atom const & a2)
{
	assert(a1.prec_ == a2.prec_);
	for ( int i = 3; i-- ; )
	if ( std::abs(a1.pos_[i] - a2.pos_[i]) >= a1.prec_ )
		return a1.pos_[i] < a2.pos_[i];
	return false;
}

} /* namespace LatticeStructure */
} /* namespace elephon */
