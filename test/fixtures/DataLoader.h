/*	This file DataLoader.h is part of elephon.
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
 *  Created on: Jun 21, 2017
 *      Author: A. Linscheid
 */

#ifndef TEST_FIXTURES_DATALOADER_H_
#define TEST_FIXTURES_DATALOADER_H_

#include "IOMethods/VASPInterface.h"
#include <string>
#include <memory>

namespace test
{
namespace fixtures
{

class DataLoader
{
public:

	std::shared_ptr<elephon::IOMethods::VASPInterface>
		create_vasp_loader(	std::string const & contentInputFile,
							std::string fileName = std::string()) const;

	elephon::LatticeStructure::UnitCell
		load_unit_cell(		std::string const & contentInputFile,
							std::string fileName = std::string()) const;

private:

	void process_fileName(std::string & fileName ) const;
};

} /* namespace fixtures */
} /* namespace test */

#endif /* TEST_FIXTURES_DATALOADER_H_ */
