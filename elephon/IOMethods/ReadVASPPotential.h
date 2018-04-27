/*	This file ReadVASPPotential.h is part of elephon.
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
 *  Created on: Apr 26, 2018
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_IOMETHODS_READVASPPOTENTIAL_H_
#define ELEPHON_IOMETHODS_READVASPPOTENTIAL_H_

#include <string>
#include <array>
#include <vector>

namespace elephon
{
namespace IOMethods
{

class ReadVASPPotential
{
public:

	void set_filepath(std::string filepath);

	std::string const & get_default_filename() const;

	template<class VTReg, class VTRad, class VTChg>
	void read_potential_file(
			std::array<int,3> & regularGridDim,
			VTReg & regularData,
			std::vector<int> & angularLMaxPerAtom,
			std::vector<std::vector<double>> & radialPointsPerAtom,
			std::vector<double> & radiusPerAtom,
			std::vector<VTRad> & radialData,
			std::vector<double> & coreChargeZ,
			std::vector<VTChg> & frozenCoreElectronChg) const;
private:

	std::string defaultFileName_ = "LOCPOT_AE";

	std::string filepath_;


};

} /* namespace IOMethods */
} /* namespace elephon */

#include "IOMethods/ReadVASPPotential.hpp"
#endif /* ELEPHON_IOMETHODS_READVASPPOTENTIAL_H_ */
