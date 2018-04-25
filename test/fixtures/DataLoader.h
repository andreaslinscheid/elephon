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
#include "ElectronicStructure/ElectronicBands.h"
#include "IOMethods/ResourceHandler.h"
#include <string>
#include <memory>
#include <boost/filesystem.hpp>

namespace elephon
{
namespace test
{
namespace fixtures
{

/**
 * Class that creates data that is used in unit tests.
 */
class DataLoader
{
public:

	std::shared_ptr<elephon::IOMethods::ResourceHandler>
		create_resource_handler(
				std::string const & contentInputFile ) const;

	std::shared_ptr<elephon::IOMethods::VASPInterface>
		create_vasp_loader(	std::string const & contentInputFile,
							std::string fileName = std::string()) const;

	elephon::LatticeStructure::UnitCell
		load_unit_cell(		std::string const & contentInputFile,
							std::string fileName = std::string()) const;

	elephon::ElectronicStructure::ElectronicBands
		create_symmetric_cosine_model(
				std::vector<int> griddims,
				std::vector<double> gridshift) const;

	elephon::LatticeStructure::Symmetry create_partial_sym() const;

	/**
	 * Create a symmetry group that features 90deg rotations
	 * about all 3 axis for a trivial lattice of length 1 in each direction.
	 *
	 * @return The symmetry group object representing this group.
	 */
	elephon::LatticeStructure::Symmetry create_4fold_symmetry() const;

	/**
	 * Get hard coded force data from vasp for the 4x4x4 primitive supercell
	 * @return	a vector set to the force data.
	 */
	std::vector<double> get_reference_force_data_vasp_Al_sc4x4x4() const;
private:

	void process_fileName(std::string & fileName ) const;
};

} /* namespace fixtures */
} /* namespace test */
} /* namespace elephon */

#endif /* TEST_FIXTURES_DATALOADER_H_ */
