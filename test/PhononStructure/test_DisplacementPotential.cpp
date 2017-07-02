/*	This file test_DisplacementPotential.cpp is part of elephon.
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
 *  Created on: Jul 1, 2017
 *      Author: A. Linscheid
 */


#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE PhononStructure
#include <boost/test/unit_test.hpp>
#include "fixtures/MockStartup.h"
#include "fixtures/DataLoader.h"
#include "PhononStructure/DisplacementPotential.h"

BOOST_AUTO_TEST_CASE( build_Al_fcc_primitive )
{
	test::fixtures::MockStartup ms;
	auto rootDir = ms.get_data_for_testing_dir() / "Al" / "vasp" / "fcc_primitive" ;
	auto phononDir = rootDir / "phonon";

	std::shared_ptr<elephon::IOMethods::VASPInterface> loader;
	std::vector<elephon::LatticeStructure::Atom> atomsUC;
	elephon::LatticeStructure::Symmetry symmetry;
	elephon::LatticeStructure::RegularGrid kgrid;
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
	elephon::LatticeStructure::RegularGrid rsGridSC;
	rsGridSC.initialize( 1e-6,
			dims,
			std::vector<double>({0.0,0.0,0.0}), //no shift
			supercell.get_symmetry(),
			supercell.get_lattice() );

	//Read the normal periodic potential
	loader->read_electronic_potential(
			rootDir.string(),
			dims,
			thisDisplPot);
	elephon::LatticeStructure::RegularGrid rsGridUC;
	rsGridUC.initialize( 1e-6,
			dims,
			std::vector<double>({0.0,0.0,0.0}), //no shift
			unitCell.get_symmetry(),
			unitCell.get_lattice() );

	elephon::PhononStructure::DisplacementPotential dvscf;
	dvscf.build( unitCell, supercell, irreducibleDispl, rsGridUC, rsGridSC,
			thisDisplPot, displPot);

	BOOST_REQUIRE( dvscf.get_num_modes() == 3 );

	BOOST_REQUIRE( dvscf.get_num_R() == 2*2*2 );
}
