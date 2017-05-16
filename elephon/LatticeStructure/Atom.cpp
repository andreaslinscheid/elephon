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

namespace elephon
{
namespace LatticeStructure
{

Atom::Atom(std::string kind, std::vector<double> pos,
		std::vector<bool> frozen)
		: kind_(std::move(kind)),pos_(std::move(pos)),frozen_(std::move(frozen))
{

};

std::string Atom::get_kind() const
{
	return kind_;
};

std::vector<double> Atom::get_position() const
{
	return pos_;
};

void Atom::set_position(std::vector<double> newPos)
{
	pos_ = std::move(newPos);
}

std::vector<bool> Atom::get_movement_fixed() const
{
	return frozen_;
};

} /* namespace LatticeStructure */
} /* namespace elephon */
