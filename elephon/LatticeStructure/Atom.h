/*	This file Atom.h is part of elephon.
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

#ifndef ELEPHON_LATTICESTRUCTURE_ATOM_H_
#define ELEPHON_LATTICESTRUCTURE_ATOM_H_

#include <string>
#include <vector>

namespace elephon
{
namespace LatticeStructure
{

class Atom
{
public:
	Atom(std::string kind, std::vector<double> pos,std::vector<bool> frozen);

	std::string get_kind() const;

	std::vector<double> get_position() const;

	void set_position(std::vector<double> newPos);

	std::vector<bool> get_movement_fixed() const;

private:

	std::string kind_;

	std::vector<double> pos_;

	std::vector<bool> frozen_;
};

} /* namespace LatticeStructure */
} /* namespace elephon */

#endif /* ELEPHON_LATTICESTRUCTURE_ATOM_H_ */
