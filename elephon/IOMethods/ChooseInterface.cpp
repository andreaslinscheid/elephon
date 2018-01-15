/*	This file ChooseInterface.cpp is part of elephon.
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
 *  Created on: Jul 24, 2017
 *      Author: A. Linscheid
 */

#include "ChooseInterface.h"
#include "IOMethods/VASPInterface.h"
#include "boost/filesystem.hpp"
#include <stdexcept>
#include <set>

namespace elephon
{
namespace IOMethods
{

std::shared_ptr<ElectronicStructureCodeInterface>
choose_interface( InputOptions const & inputOPts )
{
	namespace fs = boost::filesystem;
	auto d = inputOPts.get_root_dir();
	auto bp = fs::path(d);

	if ( not fs::exists(bp) )
		throw std::runtime_error("Directory "+d+" is not valid");

	std::set<std::string> fileSet;
	fs::recursive_directory_iterator end;

	for (fs::recursive_directory_iterator i(bp); i != end; ++i)
	{
		const fs::path cp = (*i);
		fileSet.insert(cp.filename().string());
	}

	std::shared_ptr<ElectronicStructureCodeInterface> result = nullptr;

	//Check if the files match VASP
	auto e = fileSet.end();
	if ( (fileSet.find("POSCAR") != e) or (fileSet.find("vasprun.xml") != e) or
		 (fileSet.find("POTCAR") != e) or (fileSet.find("INCAR")   != e) )
		result = std::make_shared< elephon::IOMethods::VASPInterface >(inputOPts);

	if ( result == nullptr )
	{
		std::cout << "Unable to match files in " << d <<" to an electronic structure code. Exiting ..." << std::endl;
		std::exit(0);
	}
	return result;
};

} /* namespace IOMethods */
} /* namespace elephon */
