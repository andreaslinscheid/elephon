/*	This file DataLoader.cpp is part of elephon.
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

#include <boost/filesystem.hpp>
#include "fixtures/DataLoader.h"
#include "fixtures/MockStartup.h"
#include "LatticeStructure/UnitCell.h"
#include "LatticeStructure/Symmetry.h"
#include "LatticeStructure/RegularGrid.h"
#include "LatticeStructure/LatticeModule.h"
#include "LatticeStructure/Atom.h"

namespace test
{
namespace fixtures
{

std::shared_ptr<elephon::IOMethods::VASPInterface>
DataLoader::create_vasp_loader(	std::string const & contentInputFile,
								std::string fileName) const
{
	this->process_fileName(fileName);
	elephon::IOMethods::InputOptions options;
	test::fixtures::MockStartup ms;
	ms.simulate_elephon_input( fileName, contentInputFile, options );

	return std::make_shared< elephon::IOMethods::VASPInterface >(options);
}

elephon::LatticeStructure::UnitCell
DataLoader::load_unit_cell(		std::string const & contentInputFile,
								std::string fileName ) const
{
	this->process_fileName(fileName);
	auto loader = this->create_vasp_loader(contentInputFile,fileName);

	elephon::IOMethods::InputOptions options;
	test::fixtures::MockStartup ms;
	ms.simulate_elephon_input( fileName, contentInputFile, options );

	std::vector<elephon::LatticeStructure::Atom> atomsUC;
	elephon::LatticeStructure::Symmetry symmetry;
	elephon::LatticeStructure::RegularGrid kgrid;
	elephon::LatticeStructure::LatticeModule  lattice;

	loader->read_cell_paramters( options.get_root_dir(), options.get_gPrec(),
			kgrid,lattice,atomsUC, symmetry);

	elephon::LatticeStructure::UnitCell uc;
	uc.initialize( atomsUC, lattice, symmetry );

	return uc;
}

void
DataLoader::process_fileName(std::string & fileName ) const
{
	if ( fileName.empty() )
			fileName = (boost::filesystem::temp_directory_path() / boost::filesystem::unique_path()).string();
}

} /* namespace fixtures */
} /* namespace test */
