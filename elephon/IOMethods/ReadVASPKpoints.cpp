/*	This file ReadVASPKpoints.cpp is part of elephon.
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
 *  Created on: Sep 27, 2017
 *      Author: A. Linscheid
 */

#include "IOMethods/ReadVASPKpoints.h"
#include <fstream>
#include <stdexcept>
#include <vector>
#include <cmath>

namespace elephon
{
namespace IOMethods
{

void
ReadVASPKpoints::read_kpoints(
		std::string filename,
		LatticeStructure::LatticeModule const & lattice)
{
	filename_ = std::move(filename);

	std::ifstream file( filename_.c_str() );
	if ( ! file.good())
		throw std::runtime_error(std::string("Failed to open KPOINTS file for reading: ") + filename_);

	std::string line;

	// comment
	std::getline( file, line);

	// switch automatic
	std::getline( file, line);
	int numberOfK = std::stoi(line);
	if ( numberOfK != 0 )
		throw std::runtime_error("Reading of other than automatic meshes not implemented");

	// switch Monkhurst mesh
	std::getline( file, line);
	if ( line.empty() )
		throw std::runtime_error(std::string("Incorrect format in file ")+filename_);

	if ( (line.front() == 'g') or (line.front() == 'G')
			or (line.front() == 'm') or (line.front() == 'M'))
	{
		bool monk = (line.front() == 'm') or (line.front() == 'M');
		std::getline( file, line);
		size_t sz;
		dim_[0] = std::stoi(line,&sz);
		auto remainder = line.substr(sz);
		dim_[1] = std::stoi(remainder,&sz);
		remainder = remainder.substr(sz);
		dim_[2] = std::stoi(remainder,&sz);
		if ( monk )
			for (int i = 0 ; i < 3 ; ++i)
				gridShift_[i] = 0.5*double((dim_[i]+1) % 2); // even grids get a shift of 0.5

		if ( std::getline( file, line) )
		{
			gridShift_[0] += std::stoi(line,&sz);
			remainder = line.substr(sz);
			gridShift_[1] += std::stoi(remainder,&sz);
			remainder = remainder.substr(sz);
			gridShift_[2] += std::stoi(remainder,&sz);
		}
	}
	if ( (line.front() == 'a') or (line.front() == 'A') )
	{
		std::getline( file, line);
		double length = std::stof(line);

		auto norm = [] (std::vector<double> const & b ) { return std::sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);};
		// shift = 0
		auto b1 = lattice.get_reci_lattice_vector(0);
		auto b2 = lattice.get_reci_lattice_vector(1);
		auto b3 = lattice.get_reci_lattice_vector(2);
		dim_[0] = std::floor(std::max(1.0, length*norm(b1)/lattice.get_alat() + 0.5));
		dim_[1] = std::floor(std::max(1.0, length*norm(b2)/lattice.get_alat() + 0.5));
		dim_[2] = std::floor(std::max(1.0, length*norm(b3)/lattice.get_alat() + 0.5));
	}

	// map back to lattice shifts [0, 1[
	for ( auto & s : gridShift_ )
		s -= std::floor(s);
}

std::vector<double> const &
ReadVASPKpoints::get_grid_shift() const
{
	return gridShift_;
}

std::vector<int> const &
ReadVASPKpoints::get_grid_dim() const
{
	return dim_;
}

} /* namespace IOMethods */
} /* namespace elephon */
