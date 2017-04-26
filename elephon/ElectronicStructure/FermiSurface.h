/*	This file FermiSurface.h is part of elephon.
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
 *  Created on: Apr 26, 2017
 *      Author: A. Linscheid
 */

#ifndef ELECTRONICSTRUCTURE_FERMISURFACE_H_
#define ELECTRONICSTRUCTURE_FERMISURFACE_H_

#include <vector>
#include <cstdlib>

namespace elephon
{
namespace ElectronicStructure
{

class FermiSurface
{
public:
	FermiSurface();

	void triangulate_Fermi_surface(
			std::vector<size_t> const& kgrid,
			std::vector<double> const& energies,
			size_t numTargetPoints);

	size_t get_npts() const;

	void get_pt(size_t i, std::vector<double> & p) const;

	void get_pt_weight(size_t i, double & pw) const;
private:

	std::vector<double> kfPoints_;

	std::vector<double> kfWeights_;
};

} /* namespace ElectronicStructure */
} /* namespace elephon */

#endif /* ELECTRONICSTRUCTURE_FERMISURFACE_H_ */
