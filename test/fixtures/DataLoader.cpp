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

elephon::LatticeStructure::Symmetry
DataLoader::create_4fold_symmetry() const
{
	elephon::LatticeStructure::Symmetry sym;
	std::vector<int> symops{  1,  0,  0,  0,  1,  0,  0,  0,  1, // No: 1
							 -1,  0,  0,  0, -1,  0,  0,  0, -1, // No: 2
							  0,  0,  1,  1,  0,  0,  0,  1,  0, // No: 3
							  0,  0, -1, -1,  0,  0,  0, -1,  0, // No: 4
							  0,  1,  0,  0,  0,  1,  1,  0,  0, // No: 5
							  0, -1,  0,  0,  0, -1, -1,  0,  0, // No: 6
							  0, -1,  0,  1,  0,  0,  0,  0,  1, // No: 7
							  0,  1,  0, -1,  0,  0,  0,  0, -1, // No: 8
							 -1,  0,  0,  0,  0,  1,  0,  1,  0, // No: 9
							  1,  0,  0,  0,  0, -1,  0, -1,  0, // No:10
							  0,  0, -1,  0,  1,  0,  1,  0,  0, // No:11
							  0,  0,  1,  0, -1,  0, -1,  0,  0, // No:12
							 -1,  0,  0,  0, -1,  0,  0,  0,  1, // No:13
							  1,  0,  0,  0,  1,  0,  0,  0, -1, // No:14
							  0,  0, -1, -1,  0,  0,  0,  1,  0, // No:15
							  0,  0,  1,  1,  0,  0,  0, -1,  0, // No:16
							  0, -1,  0,  0,  0, -1,  1,  0,  0, // No:17
							  0,  1,  0,  0,  0,  1, -1,  0,  0, // No:18
							  0,  1,  0, -1,  0,  0,  0,  0,  1, // No:19
							  0, -1,  0,  1,  0,  0,  0,  0, -1, // No:20
							  1,  0,  0,  0,  0, -1,  0,  1,  0, // No:21
							 -1,  0,  0,  0,  0,  1,  0, -1,  0, // No:22
							  0,  0,  1,  0, -1,  0,  1,  0,  0, // No:23
							  0,  0, -1,  0,  1,  0, -1,  0,  0, // No:24
							  0,  0,  1, -1,  0,  0,  0, -1,  0, // No:25
							  0,  0, -1,  1,  0,  0,  0,  1,  0, // No:26
							  0,  1,  0,  0,  0, -1, -1,  0,  0, // No:27
							  0, -1,  0,  0,  0,  1,  1,  0,  0, // No:28
							  1,  0,  0,  0, -1,  0,  0,  0, -1, // No:29
							 -1,  0,  0,  0,  1,  0,  0,  0,  1, // No:30
							  0,  0,  1,  0,  1,  0, -1,  0,  0, // No:31
							  0,  0, -1,  0, -1,  0,  1,  0,  0, // No:32
							  0,  1,  0,  1,  0,  0,  0,  0, -1, // No:33
							  0, -1,  0, -1,  0,  0,  0,  0,  1, // No:34
							  1,  0,  0,  0,  0,  1,  0, -1,  0, // No:35
							 -1,  0,  0,  0,  0, -1,  0,  1,  0, // No:36
							  0, -1,  0,  0,  0,  1, -1,  0,  0, // No:37
							  0,  1,  0,  0,  0, -1,  1,  0,  0, // No:38
							 -1,  0,  0,  0,  1,  0,  0,  0, -1, // No:39
							  1,  0,  0,  0, -1,  0,  0,  0,  1, // No:40
							  0,  0, -1,  1,  0,  0,  0, -1,  0, // No:41
							  0,  0,  1, -1,  0,  0,  0,  1,  0, // No:42
							 -1,  0,  0,  0,  0, -1,  0, -1,  0, // No:43
							  1,  0,  0,  0,  0,  1,  0,  1,  0, // No:44
							  0, -1,  0, -1,  0,  0,  0,  0, -1, // No:45
							  0,  1,  0,  1,  0,  0,  0,  0,  1, // No:46
							  0,  0, -1,  0, -1,  0, -1,  0,  0, // No:47
							  0,  0,  1,  0,  1,  0,  1,  0,  0, // No:48
	};
	std::vector<double> frac(symops.size()/3, 0.0);
	sym.initialize(1e-6, symops, frac, elephon::LatticeStructure::LatticeModule(), true);
	return sym;
}

std::vector<double>
DataLoader::get_reference_force_data_vasp_Al_sc4x4x4() const
{
	std::vector<double> result{
	        0.03393673,  0.00000000, -0.03393673,
	       -0.00956984,  0.00000000,  0.00956984,
	       -0.00065475,  0.00000000,  0.00065475,
	       -0.00915954,  0.00000000,  0.00915954,
	        0.00115708,  0.00480979,  0.00462055,
	        0.00007651, -0.00009804, -0.00011443,
	        0.00012188,  0.00010106, -0.00008150,
	       -0.00456729, -0.00472585, -0.00113419,
	        0.00003514,  0.00032140,  0.00033309,
	       -0.00002138, -0.00000125,  0.00001967,
	       -0.00033309, -0.00032140, -0.00003514,
	       -0.00001967,  0.00000125,  0.00002138,
	        0.00113419,  0.00472585,  0.00456729,
	       -0.00462055, -0.00480979, -0.00115708,
	        0.00011443,  0.00009804, -0.00007651,
	        0.00008150, -0.00010106, -0.00012188,
	       -0.00462055,  0.00480979, -0.00115708,
	        0.00011443, -0.00009804, -0.00007651,
	        0.00008150,  0.00010106, -0.00012188,
	        0.00113419, -0.00472585,  0.00456729,
	        0.00000155, -0.00017239, -0.00000155,
	        0.00009213, -0.00000246, -0.00009213,
	        0.00000538,  0.00017798, -0.00000538,
	        0.00033735, -0.00000224, -0.00033735,
	       -0.00009726, -0.00007627,  0.00005406,
	       -0.00005989,  0.00007483,  0.00009497,
	       -0.00020494, -0.00005849, -0.00011421,
	        0.00011510,  0.00005660,  0.00020308,
	        0.00019917,  0.00000000, -0.00014828,
	       -0.00059292,  0.00000000, -0.00033925,
	       -0.00026259,  0.00000000,  0.00026048,
	        0.00033542,  0.00000000,  0.00059063,
	       -0.00033309,  0.00032140, -0.00003514,
	       -0.00001967, -0.00000125,  0.00002138,
	        0.00003514, -0.00032140,  0.00033309,
	       -0.00002138,  0.00000125,  0.00001967,
	       -0.00005406, -0.00007627,  0.00009726,
	       -0.00009497,  0.00007483,  0.00005989,
	        0.00011421, -0.00005849,  0.00020494,
	       -0.00020308,  0.00005660, -0.00011510,
	       -0.00001162,  0.00000000,  0.00001162,
	       -0.00014434,  0.00000000,  0.00014434,
	        0.00021596,  0.00000000, -0.00021596,
	       -0.00014862,  0.00000000,  0.00014862,
	       -0.00005989, -0.00007483,  0.00009497,
	       -0.00020494,  0.00005849, -0.00011421,
	        0.00011510, -0.00005660,  0.00020308,
	       -0.00009726,  0.00007627,  0.00005406,
	       -0.00456729,  0.00472585, -0.00113419,
	        0.00115708, -0.00480979,  0.00462055,
	        0.00007651,  0.00009804, -0.00011443,
	        0.00012188, -0.00010106, -0.00008150,
	        0.00014828,  0.00000000, -0.00019917,
	        0.00033925,  0.00000000,  0.00059292,
	       -0.00026048,  0.00000000,  0.00026259,
	       -0.00059063,  0.00000000, -0.00033542,
	       -0.00009497, -0.00007483,  0.00005989,
	        0.00011421,  0.00005849,  0.00020494,
	       -0.00020308, -0.00005660, -0.00011510,
	       -0.00005406,  0.00007627,  0.00009726,
	        0.00000538, -0.00017798, -0.00000538,
	        0.00033735,  0.00000224, -0.00033735,
	        0.00000155,  0.00017239, -0.00000155,
	        0.00009213,  0.00000246, -0.00009213,
	};
	return result;
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
