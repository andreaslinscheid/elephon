/*	This file ReadVASPPoscar.h is part of elephon.
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
 *  Created on: May 14, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_IOMETHODS_READVASPPOSCAR_H_
#define ELEPHON_IOMETHODS_READVASPPOSCAR_H_

#include <string>
#include <vector>

namespace elephon
{
namespace IOMethods
{

/**
 *
 */
class ReadVASPPoscar
{
public:

	typedef class AtomPos_
	{
	public:
		AtomPos_(std::string kind, std::vector<double> pos,
				std::vector<bool> frozen)
				: kind_(std::move(kind)),pos_(std::move(pos)),frozen_(std::move(frozen)) {	};

		std::string get_kind() const {	return kind_; };

		std::vector<double> get_position() const {	return pos_; };

		std::vector<bool> get_movement_fixed() const {	return frozen_; };

	private:
		std::string kind_;
		std::vector<double> pos_;
		std::vector<bool> frozen_;
	} AtomPos;

	void read_file( std::string filename );

	std::vector<double> get_lattice_matrix() const;

	std::vector<AtomPos> get_atoms_list() const;

private:

	std::vector<double> latticeMatrix_;

	std::vector<AtomPos> atoms_;
};

} /* namespace IOMethods */
} /* namespace elephon */

#endif /* ELEPHON_IOMETHODS_READVASPPOSCAR_H_ */
