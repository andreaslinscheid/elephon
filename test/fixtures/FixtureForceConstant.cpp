/*	This file FixtureForceConstant.cpp is part of elephon.
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

#include "fixtures/FixtureForceConstant.h"
#include "fixtures/DataLoader.h"
#include "fixtures/MockStartup.h"
#include <boost/filesystem.hpp>

namespace test
{
namespace fixtures
{

elephon::PhononStructure::ForceConstantMatrix
FixtureForceConstant::compute_fc_for_Al_gamma()
{
	elephon::PhononStructure::ForceConstantMatrix fc;
	test::fixtures::MockStartup ms;
	auto rootDir = ms.get_data_for_testing_dir() / "Al" / "vasp" / "fcc_primitive" ;
	auto phononDir = rootDir / "phonon";

	std::string content = std::string()+
			"scell=2 2 2\n"
			"root_dir="+rootDir.string()+"\n"
			"elphd="+phononDir.string()+"\n"
			"";

	DataLoader dl;
	auto uc = dl.load_unit_cell(content);
	auto loader = dl.create_vasp_loader(content);

	//Here, we regenerate the displacement
	std::vector<elephon::LatticeStructure::AtomDisplacement> irreducibleDispl;
	uc.generate_displacements(0.01,
			true,
			irreducibleDispl);

	//Here, we read the forces from the vasp output
	int nIrdDispl = int(irreducibleDispl.size());
	std::vector<std::vector<double>> forces( nIrdDispl );
	std::vector<double> thisForces;
	for ( int idispl = 0 ; idispl < nIrdDispl; ++idispl )
	{
		loader->read_forces(
				(phononDir / (std::string("displ_")+std::to_string(idispl))).string(),
				thisForces);
		forces[idispl] = std::move(thisForces);
	}

	auto sc = uc.build_supercell( 2 , 2 , 2 );
	fc.build(uc, sc, irreducibleDispl, forces );
	return fc;
}

} /* namespace fixtures */
} /* namespace test */
