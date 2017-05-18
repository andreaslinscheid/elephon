/*	This file ElectronicStructureCodeInterface.cpp is part of elephon.
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
 *  Created on: May 17, 2017
 *      Author: A. Linscheid
 */

#include "ElectronicStructureCodeInterface.h"
#include <boost/filesystem.hpp>

namespace elephon
{
namespace IOMethods
{

ElectronicStructureCodeInterface::ElectronicStructureCodeInterface(
		IOMethods::InputOptions inputOPts ) : inputOPts_(inputOPts)
{
}

ElectronicStructureCodeInterface::~ElectronicStructureCodeInterface()
{
	// TODO this class could take care of keeping track of file
}

std::vector<std::string>
ElectronicStructureCodeInterface::gen_input_file_list(std::string directory) const
{
	boost::filesystem::path dir(directory);
	auto result = this->list_all_input_files();
	for ( auto &f : result )
		f = (dir / boost::filesystem::path(f)).string();
	return result;
}

} /* namespace IOMethods */
} /* namespace elephon */
