/*	This file PhononGrid.cpp is part of elephon.
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
 *  Created on: Oct 2, 2017
 *      Author: A. Linscheid
 */

#include "PhononStructure/PhononGrid.h"
#include <chrono>
#include <fstream>

namespace elephon
{
namespace PhononStructure
{

int
PhononGrid::num_modes() const
{
	return this->get_nData_gpt();
}

void
PhononGrid::write_phonon_dos_file(
		std::string const & filename,
		std::vector<double> const & frequencies,
		std::shared_ptr<const LatticeStructure::TetrahedraGrid> tetra,
		std::shared_ptr<const Phonon> ph) const
{
	std::ofstream file( filename.c_str() );
	if ( not file.good() )
		throw std::runtime_error(std::string("Problem opening file ")+filename+" for writing the phonon DOS data");

	std::vector<double> phdos;
	if ( ! tetra )
	{
		if ( ! ph )
			this->compute_DOS(frequencies, phdos);
		else
			this->compute_DOS_wan(*ph, frequencies, phdos);
	}
	else
		this->compute_DOS_tetra(tetra, frequencies, phdos);

	auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
	file << "# phonon DOS in units of 1/THz with frequency in units of THz. Date is " << std::ctime(&now);

	for ( int iw = 0 ; iw < frequencies.size() ; ++iw)
		file << frequencies[iw] << '\t' << phdos[iw] << '\n';
	file.close();
}

} /* namespace PhononStructure */
} /* namespace elephon */
