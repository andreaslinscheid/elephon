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
#include "fixtures/FixtureForceConstant.h"
#include "PhononStructure/DisplacementPotential.h"
#include "PhononStructure/Phonon.h"

void build_displ_pot_Al_fcc_primitive_vasp(
		elephon::PhononStructure::DisplacementPotential & dvscf)
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

	dvscf.build( unitCell, supercell, irreducibleDispl, rsGridUC, rsGridSC,
			thisDisplPot, displPot);
}

BOOST_AUTO_TEST_CASE( build_Al_fcc_primitive )
{
	elephon::PhononStructure::DisplacementPotential dvscf;
	build_displ_pot_Al_fcc_primitive_vasp(dvscf);

	BOOST_REQUIRE( dvscf.get_num_modes() == 3 );

	BOOST_REQUIRE( dvscf.get_num_R() == 2*2*2 );

	//TODO invent further checks that this is correct ...
}

BOOST_AUTO_TEST_CASE( plot_dvscf_q_Al_fcc_primitive )
{
	test::fixtures::MockStartup ms;
	auto rootDir = ms.get_data_for_testing_dir() / "Al" / "vasp" / "fcc_primitive" ;
	elephon::PhononStructure::DisplacementPotential dvscf;
	build_displ_pot_Al_fcc_primitive_vasp(dvscf);

	//Write the real space variant
	dvscf.write_dvscf(0,0,(rootDir / "dvscf.dat").string());

	//Write the q displacement variant
	std::vector<double> qVect{ 0.0,0.0,0.0 , 0.25,0.0,0.0, 0.5,0.0,0.0 };
	std::vector<int> modes{0,1};
	test::fixtures::FixtureForceConstant ffc;
	auto fc = ffc.compute_fc_for_Al_gamma();
	std::vector<double> masses = {26.9815385};

	elephon::PhononStructure::Phonon ph;
	ph.initialize( fc, masses );
	std::vector<double> w;
	std::vector< std::complex<double> > dynMat;
	ph.compute_at_q( qVect, w, dynMat );

	dvscf.write_dvscf_q(qVect,modes,dynMat,masses,(rootDir / "dvscf_q.dat").string());

	//At this point we can perform tests on the files and its content.

	//Outcomment the following for manual inspection of the files generated
	BOOST_REQUIRE( boost::filesystem::is_regular_file(rootDir / "dvscf.dat") );
	boost::filesystem::remove( rootDir / "dvscf.dat" );

	for ( auto mu : modes )
	{
		for ( int iq = 0 ; iq < qVect.size()/3; ++iq)
		{
			std::string filename = std::string("dvscf_q_")+std::to_string(iq)+"_"+std::to_string(mu)+".dat" ;
			BOOST_REQUIRE( boost::filesystem::is_regular_file( rootDir / filename ) );
			boost::filesystem::remove( rootDir / filename );
		}
	}
}
