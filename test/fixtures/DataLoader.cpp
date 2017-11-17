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
#include <LatticeStructure/RegularSymmetricGrid.h>
#include "fixtures/DataLoader.h"
#include "fixtures/MockStartup.h"
#include "LatticeStructure/UnitCell.h"
#include "LatticeStructure/Symmetry.h"
#include "LatticeStructure/LatticeModule.h"
#include "LatticeStructure/Atom.h"

namespace elephon
{
namespace test
{
namespace fixtures
{


std::shared_ptr<elephon::IOMethods::ResourceHandler>
DataLoader::create_resource_handler(
			std::string const & contentInputFile ) const
{
	auto loader = this->create_vasp_loader(contentInputFile);
	auto res = elephon::IOMethods::ResourceHandler(loader);
	return std::make_shared<elephon::IOMethods::ResourceHandler>(std::move(res));
}

std::shared_ptr<elephon::IOMethods::VASPInterface>
DataLoader::create_vasp_loader(	std::string const & contentInputFile,
								std::string fileName) const
{
	this->process_fileName(fileName);
	elephon::IOMethods::InputOptions options;
	MockStartup ms;
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
	elephon::LatticeStructure::RegularSymmetricGrid kgrid;
	elephon::LatticeStructure::LatticeModule  lattice;

	loader->read_cell_paramters( options.get_root_dir(), options.get_gPrec(),
			kgrid,lattice,atomsUC, symmetry);

	elephon::LatticeStructure::UnitCell uc;
	uc.initialize( atomsUC, lattice, symmetry );

	return uc;
}

elephon::ElectronicStructure::ElectronicBands
DataLoader::create_symmetric_cosine_model(
		std::vector<int> griddims,
		std::vector<double> gridshift) const
{
	int nB = 2;

	Auxillary::alignedvector::DV data(nB*griddims[2]*griddims[1]*griddims[0]);
	for ( int iz = 0 ; iz < griddims[2]; ++iz )
		for ( int iy = 0 ; iy < griddims[1]; ++iy )
			for ( int ix = 0 ; ix < griddims[0]; ++ix )
			{
				int cnsq = ix + griddims[0]*(iy + griddims[1]*iz);
				data[cnsq*nB+0] = std::cos((2*M_PI/griddims[0])*(ix+gridshift[0]))
								 +std::cos((2*M_PI/griddims[1])*(iy+gridshift[1]));

				data[cnsq*nB+1] = std::cos((2*M_PI/griddims[2])*(iz+gridshift[2]));
			}

	elephon::ElectronicStructure::ElectronicBands bands;
	elephon::LatticeStructure::Symmetry sym = this->create_partial_sym();
	sym.set_reciprocal_space_sym();
	elephon::LatticeStructure::RegularSymmetricGrid grid;
	grid.initialize(
			griddims,
			1e-6,
			gridshift,
			sym,
			elephon::LatticeStructure::LatticeModule());
	bands.initialize(nB, 0.0, data, grid);
	return bands;
}

elephon::LatticeStructure::Symmetry
DataLoader::create_partial_sym() const
{
	elephon::LatticeStructure::Symmetry sym;
	std::vector<int> symops{ 	1,  0,  0,  0,  1,  0,  0,  0,  1, // identity
							   -1,  0,  0,  0, -1,  0,  0,  0, -1, // inversion
								0,  1,  0, -1,  0,  0,  0,  0,  1, // rotation 90deg
								0, -1,  0,  1,  0,  0,  0,  0,  1, // rotation -90deg
								0, -1,  0,  1,  0,  0,  0,  0, -1, // rotation 90deg + inv
								0,  1,  0, -1,  0,  0,  0,  0, -1, // rotation 90deg + inv
								1,  0,  0,  0,  1,  0,  0,  0, -1,
							   -1,  0,  0,  0, -1,  0,  0,  0,  1,};// rotation 180 deg
	std::vector<double> frac(symops.size()/3, 0.0);
	sym.initialize(1e-6, symops, frac, elephon::LatticeStructure::LatticeModule(), true);
	return sym;
}

void
DataLoader::process_fileName(std::string & fileName ) const
{
	if ( fileName.empty() )
			fileName = (boost::filesystem::temp_directory_path() / boost::filesystem::unique_path()).string();
}

} /* namespace fixtures */
} /* namespace test */
} /* namespace elephon */
