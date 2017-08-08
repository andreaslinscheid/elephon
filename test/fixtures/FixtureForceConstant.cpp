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
	auto rootDir = ms.get_data_for_testing_dir() / "Al" / "vasp" / "fcc_primitive" / "phonon_run" ;
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

elephon::PhononStructure::DisplacementPotential
FixtureForceConstant::build_displ_pot_Al_fcc_primitive_vasp_sc2x2x2( )
{
	test::fixtures::MockStartup ms;
	auto rootDir = ms.get_data_for_testing_dir() / "Al" / "vasp" / "fcc_primitive" /  "phonon_run";
	auto phononDir = rootDir / "phonon";

	std::shared_ptr<elephon::IOMethods::VASPInterface> loader;
	std::vector<elephon::LatticeStructure::Atom> atomsUC;
	elephon::LatticeStructure::Symmetry symmetry;
	elephon::LatticeStructure::RegularSymmetricGrid kgrid;
	elephon::LatticeStructure::LatticeModule  lattice;
	//here we create the test input file
	std::string content = std::string()+
			"scell=2 2 2\n"
			"root_dir="+rootDir.string()+"\n"
			"elphd="+phononDir.string()+"\n"
			"";

	test::fixtures::DataLoader dl;
	loader = dl.create_vasp_loader( content );

	loader->read_cell_paramters( rootDir.string(),
			1e-6,kgrid, lattice, atomsUC, symmetry);

	elephon::LatticeStructure::UnitCell unitCell;
	unitCell.initialize(atomsUC,lattice, symmetry);

	//here we build the supercell that was used to generate the test data.
	auto supercell = unitCell.build_supercell( 2, 2, 2 );

	//Here, we regenerate the displacement
	std::vector<elephon::LatticeStructure::AtomDisplacement> irreducibleDispl;
	unitCell.generate_displacements(0.01,
			true,
			irreducibleDispl);

	//Here, we read the potential from the vasp output
	int nIrdDispl = int(irreducibleDispl.size());
	std::vector<std::vector<double>> displPot( nIrdDispl );
	std::vector<double> thisDisplPot;
	std::vector<int> dims;
	for ( int idispl = 0 ; idispl < nIrdDispl; ++idispl )
	{
		loader->read_electronic_potential(
				(phononDir / (std::string("displ_")+std::to_string(idispl))).string(),
				dims,
				thisDisplPot);

		displPot[idispl] = std::move(thisDisplPot);
	}
	elephon::LatticeStructure::RegularSymmetricGrid rsGridSC;
	rsGridSC.initialize( dims );

	//Read the normal periodic potential
	loader->read_electronic_potential(
			rootDir.string(),
			dims,
			thisDisplPot);
	elephon::LatticeStructure::RegularSymmetricGrid rsGridUC;
	rsGridUC.initialize( dims );

	elephon::PhononStructure::DisplacementPotential dvscf;
	dvscf.build( unitCell, supercell, irreducibleDispl, rsGridUC, rsGridSC,
			thisDisplPot, displPot);
	return dvscf;
}

} /* namespace fixtures */
} /* namespace test */
