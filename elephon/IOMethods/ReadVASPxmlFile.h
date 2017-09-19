/*	This file ReadVASPxmlFile.h is part of elephon.
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
 *  Created on: May 31, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_IOMETHODS_READVASPXMLFILE_H_
#define ELEPHON_IOMETHODS_READVASPXMLFILE_H_

#include <vector>
#include <string>

namespace elephon
{
namespace IOMethods
{

class ReadVASPxmlFile
{
public:

	bool is_parsed() const;

	void parse_file( std::string filename );

	std::vector<double> const & get_forces() const;

	std::vector<double> const & get_k_points() const;

	double get_Fermi_energy() const;

	std::vector<double> const & get_energies() const;

	int get_nBnd() const;

	int get_nkp() const;
private:

	int nBnd_ = 0;

	double eFermi_ = 0;

	std::string filename_;

	std::vector<double> forces_;

	std::vector<double> kpoints_;

	std::vector<double> energies_;

	std::vector<double> latticeMat_;

	std::vector<std::vector<double>> atomicPos_;
};

} /* namespace IOMethods */
} /* namespace elephon */

#endif /* ELEPHON_IOMETHODS_READVASPXMLFILE_H_ */
