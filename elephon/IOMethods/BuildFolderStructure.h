/*	This file BuildFolderStructure.h is part of elephon.
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

#ifndef ELEPHON_IOMETHODS_BUILDFOLDERSTRUCTURE_H_
#define ELEPHON_IOMETHODS_BUILDFOLDERSTRUCTURE_H_

#include "IOMethods/InputOptions.h"
#include "IOMethods/ElectronicStructureCodeInterface.h"
#include "IOMethods/ResourceHandler.h"
#include "LatticeStructure/AtomDisplacement.h"
#include "LatticeStructure/UnitCell.h"
#include <memory>
#include <string>
#include <boost/filesystem.hpp>

namespace elephon
{
namespace IOMethods
{

class BuildFolderStructure
{
public:

	void build(std::shared_ptr<ResourceHandler> resources) const;

	void build(
			IOMethods::InputOptions const & input,
			LatticeStructure::UnitCell const & unitcell,
			ElectronicStructureCodeInterface & interface ) const;

	bool check_is_build( std::string rootFolder ) const;

private:
};

} /* namespace IOMethods */
} /* namespace elephon */

#endif /* ELEPHON_IOMETHODS_BUILDFOLDERSTRUCTURE_H_ */
