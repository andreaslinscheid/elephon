/*	This file test_ForceConstantMatrix.cpp is part of elephon.
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
 *  Created on: Jun 1, 2017
 *      Author: A. Linscheid
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Input_test
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "PhononStructure/ForceConstantMatrix.h"
#include "IOMethods/Input.h"
#include "IOMethods/VASPInterface.h"
#include <assert.h>
#include <vector>
#include <fstream>
#include <string>

void load_data(
		boost::filesystem::path rootDir ,
		std::shared_ptr<elephon::IOMethods::VASPInterface> & loader,
		std::vector<elephon::LatticeStructure::Atom> & atomsUC,
		elephon::LatticeStructure::Symmetry & symmetry,
		elephon::LatticeStructure::RegularGrid & kgrid,
		elephon::LatticeStructure::LatticeModule & lattice)
{
	//here we create the test input file
	std::string content = std::string()+
			"scell=1 1 1\n"
			"root_dir="+rootDir.string()+"\n"
			"elphd="+rootDir.string()+"\n"
			"";
	std::ofstream file( (rootDir / "test_elephon_input.dat").c_str() );
	file << content;
	file.close();

	//here we read the input file via elephons input mechanism
	std::string progn("program name");
	std::string fileName = (rootDir / "test_elephon_input.dat").string();
	char * prog = new char [progn.size()+1];
	std::copy(progn.begin(),progn.end(),prog);
	prog[progn.size()] = '\0';
	char * arg =  new char [fileName.size()+1];
	std::copy(fileName.begin(),fileName.end(),arg);
	arg[fileName.size()] = '\0';
	char *argv[] = {prog, arg, NULL};
	int argc = sizeof(argv) / sizeof(char*) - 1;
	elephon::IOMethods::Input input(argc,argv);
	elephon::IOMethods::InputOptions options = input.get_opts();
	delete [] prog;
	delete [] arg;

	//introduce the VASP data loader
	loader = std::make_shared< elephon::IOMethods::VASPInterface >(input.get_opts());

	loader->read_cell_paramters( rootDir.string(),
			1e-6,kgrid, lattice, atomsUC, symmetry);
}

/**
 * Instead of loading the symmetry operation from the data in rootDir
 * set it explicitly with symmetry.
 */
void run_test(boost::filesystem::path rootDir ,
		std::shared_ptr<elephon::IOMethods::VASPInterface> loader,
		std::vector<elephon::LatticeStructure::Atom> const & atomsUC,
		elephon::LatticeStructure::Symmetry const & symmetry,
		elephon::LatticeStructure::RegularGrid  const& kgrid,
		elephon::LatticeStructure::LatticeModule  const& lattice,
		std::vector<double> const & ref_forces_displ_0)
{
	elephon::LatticeStructure::UnitCell unitCell;
	unitCell.initialize(atomsUC,lattice, symmetry);

	//here we build the supercell that was used to generate the test data.
	auto supercell = unitCell.build_supercell( 1, 1, 1 );

	elephon::LatticeStructure::UnitCell reducedSymmetryUC;
	reducedSymmetryUC.initialize( unitCell.get_atoms_list(), unitCell.get_lattice(), supercell.get_symmetry());

	//Here, we regenerate the displacement
	std::vector<elephon::LatticeStructure::AtomDisplacement> irreducibleDispl;
	reducedSymmetryUC.generate_displacements(0.01,
			true,
			irreducibleDispl);

	//Here, we read the forces from the vasp output
	int nIrdDispl = int(irreducibleDispl.size());
	std::vector<std::vector<double>> forces( nIrdDispl+1 );
	std::vector<double> thisForces;
	for ( int idispl = 0 ; idispl < nIrdDispl; ++idispl )
	{
		loader->read_forces(
				(rootDir / (std::string("displ_")+std::to_string(idispl))).string(),
				thisForces);

		forces[idispl] = std::move(thisForces);
	}
	//The last element is set of residual forces in the clean supercell
	loader->read_forces(
					(rootDir / "scell").string(),
					thisForces);
	forces[nIrdDispl] = std::move(thisForces);

	elephon::PhononStructure::ForceConstantMatrix phi;
	phi.build( unitCell, supercell, irreducibleDispl, forces);

	BOOST_REQUIRE( phi.get_num_modes() == 4*3 );

	BOOST_REQUIRE( phi.get_num_R() == 2*2*2 );

	//Check the R = 0 first column
	//We can't be too greedy here - the forces are not more than 2% accurate
	BOOST_CHECK_CLOSE( phi(0, 0, 0,  0, 0) , ref_forces_displ_0[ 0] , 2 );
	BOOST_CHECK_CLOSE( phi(0, 0, 0,  1, 0) , ref_forces_displ_0[ 1] , 2 );
	BOOST_CHECK_CLOSE( phi(0, 0, 0,  2, 0) , ref_forces_displ_0[ 2] , 2 );
	BOOST_CHECK_CLOSE( phi(0, 0, 0,  3, 0) , ref_forces_displ_0[ 3] , 2 );
	BOOST_CHECK_CLOSE( phi(0, 0, 0,  4, 0) , ref_forces_displ_0[ 4] , 2 );
	BOOST_CHECK_CLOSE( phi(0, 0, 0,  5, 0) , ref_forces_displ_0[ 5] , 2 );
	BOOST_CHECK_CLOSE( phi(0, 0, 0,  6, 0) , ref_forces_displ_0[ 6] , 2 );
	BOOST_CHECK_CLOSE( phi(0, 0, 0,  7, 0) , ref_forces_displ_0[ 7] , 2 );
	BOOST_CHECK_CLOSE( phi(0, 0, 0,  8, 0) , ref_forces_displ_0[ 8] , 2 );
	BOOST_CHECK_CLOSE( phi(0, 0, 0,  9, 0) , ref_forces_displ_0[ 9] , 2 );
	BOOST_CHECK_CLOSE( phi(0, 0, 0, 10, 0) , ref_forces_displ_0[10] , 2 );
	BOOST_CHECK_CLOSE( phi(0, 0, 0, 11, 0) , ref_forces_displ_0[11] , 2 );

	//it has to equal the row too
	BOOST_CHECK_CLOSE( phi(0, 0, 0, 0,  0) , ref_forces_displ_0[ 0] , 2 );
	BOOST_CHECK_CLOSE( phi(0, 0, 0, 0,  1) , ref_forces_displ_0[ 1] , 2 );
	BOOST_CHECK_CLOSE( phi(0, 0, 0, 0,  2) , ref_forces_displ_0[ 2] , 2 );
	BOOST_CHECK_CLOSE( phi(0, 0, 0, 0,  3) , ref_forces_displ_0[ 3] , 2 );
	BOOST_CHECK_CLOSE( phi(0, 0, 0, 0,  4) , ref_forces_displ_0[ 4] , 2 );
	BOOST_CHECK_CLOSE( phi(0, 0, 0, 0,  5) , ref_forces_displ_0[ 5] , 2 );
	BOOST_CHECK_CLOSE( phi(0, 0, 0, 0,  6) , ref_forces_displ_0[ 6] , 2 );
	BOOST_CHECK_CLOSE( phi(0, 0, 0, 0,  7) , ref_forces_displ_0[ 7] , 2 );
	BOOST_CHECK_CLOSE( phi(0, 0, 0, 0,  8) , ref_forces_displ_0[ 8] , 2 );
	BOOST_CHECK_CLOSE( phi(0, 0, 0, 0,  9) , ref_forces_displ_0[ 9] , 2 );
	BOOST_CHECK_CLOSE( phi(0, 0, 0, 0, 10) , ref_forces_displ_0[10] , 2 );
	BOOST_CHECK_CLOSE( phi(0, 0, 0, 0, 11) , ref_forces_displ_0[11] , 2 );

	//While the data was generated without symmetry, it must still be approximately there
	elephon::LatticeStructure::UnitCell symmetricUC;
	symmetricUC.initialize( atomsUC, lattice, symmetry);
	std::vector<std::vector<int> > rotMap;
 	symmetricUC.generate_rotation_maps( rotMap );
	std::vector<double> localMatrix(9);
	std::vector<double> rotLocalMatrix(9);
	for ( int ia1 = 0 ; ia1 < symmetricUC.get_atoms_list().size() ; ++ia1)
		for ( int ia2 = 0 ; ia2 < symmetricUC.get_atoms_list().size() ; ++ia2)
		{
			//copy the local matrix
			for ( int i = 0 ; i < 3 ; ++i)
				for ( int j = 0 ; j < 3 ; ++j)
					localMatrix[i*3+j] = phi(0, 0, 0, ia1*3+i, ia2*3+j);

			//compare the rotated local copy with the value from the ForceConstantMatrix object
			for ( int isym = 0 ; isym < symmetricUC.get_symmetry().get_num_symmetries(); ++isym)
			{
				rotLocalMatrix = localMatrix;
				symmetricUC.get_symmetry().rotate_matrix_cartesian(isym,rotLocalMatrix.begin(),rotLocalMatrix.end());

				int ia1r = rotMap[isym][ia1];
				int ia2r = rotMap[isym][ia2];

				for ( int i = 0 ; i < 3 ; ++i)
					for ( int j = 0 ; j < 3 ; ++j)
						BOOST_CHECK_CLOSE( phi(0, 0, 0, ia1r*3+i, ia2r*3+j) , rotLocalMatrix[i*3+j] , 2 );
			}
		}
}

BOOST_AUTO_TEST_CASE( Assemble_ForceConstantMatrix_no_symmetry )
{
	boost::filesystem::path rootDir = boost::filesystem::path(__FILE__).parent_path()
			/ "../data_for_testing/Al/vasp/build_mat_fc_no_sym/";

	//In this example we don't include symmetries
	elephon::LatticeStructure::Symmetry identityGroup;

	//reference data - pulled from the file - for irreducible displacement 0
	// which displaces atom 0 in direction x and is thus element mu2 = 0
	//Use it to compute the first column of phi explicitly
	std::vector<double> ref_forces_displ_0 =	//	Positions in unperturbed cell (mapped to the 1 UC):
	{-0.02181663 ,  0.00000000 , -0.00000000, // 0.00000000  0.00000000  0.00000000  R= 0 0 0 Atom : 0
	  0.02968642 , -0.00000000 ,  0.00000000, // 0.00000000 -0.25000000 -0.50000000  R= 0 0 0 Atom : 1
	 -0.00393490 , -0.00000000 ,  0.00000000, //-0.25000000  0.00000000 -0.50000000  R= 0 0 0 Atom : 2
	 -0.00393490 , -0.00000000 ,  0.00000000};//-0.25000000 -0.25000000  0.00000000  R= 0 0 0 Atom : 3

	double magn = -0.00123800*8.07786000;

	for ( auto &f : ref_forces_displ_0)
		f /= magn;

	std::shared_ptr<elephon::IOMethods::VASPInterface> loader;
	std::vector<elephon::LatticeStructure::Atom> atomsUC;
	elephon::LatticeStructure::Symmetry symmetry;
	elephon::LatticeStructure::RegularGrid kgrid;
	elephon::LatticeStructure::LatticeModule  lattice;
	load_data( rootDir, loader, atomsUC, symmetry, kgrid, lattice);

	run_test( rootDir, loader, atomsUC, identityGroup, kgrid, lattice, ref_forces_displ_0 );
}

BOOST_AUTO_TEST_CASE( Assemble_ForceConstantMatrix_partial_symmetry )
{
	boost::filesystem::path rootDir = boost::filesystem::path(__FILE__).parent_path()
			/ "../data_for_testing/Al/vasp/build_mat_fc_partial_sym/";

	//reference data - pulled from the file - for irreducible displacement 0
	// which displaces atom 0 in direction x and is thus element mu2 = 0
	//Use it to compute the first column of phi explicitly
	std::vector<double> ref_forces_displ_0 =	//	Positions in unperturbed cell (mapped to the 1 UC):
	{-0.02181663 ,  0.00000000 , -0.00000000, // 0.00000000  0.00000000  0.00000000  R= 0 0 0 Atom : 0
	  0.02968642 , -0.00000000 ,  0.00000000, // 0.00000000 -0.25000000 -0.50000000  R= 0 0 0 Atom : 1
	 -0.00393490 , -0.00000000 ,  0.00000000, //-0.25000000  0.00000000 -0.50000000  R= 0 0 0 Atom : 2
	 -0.00393490 , -0.00000000 ,  0.00000000};//-0.25000000 -0.25000000  0.00000000  R= 0 0 0 Atom : 3

	double magn = -0.00123800*8.07786000;

	for ( auto &f : ref_forces_displ_0)
		f /= magn;

	std::shared_ptr<elephon::IOMethods::VASPInterface> loader;
	std::vector<elephon::LatticeStructure::Atom> atomsUC;
	elephon::LatticeStructure::Symmetry symmetry;
	elephon::LatticeStructure::RegularGrid kgrid;
	elephon::LatticeStructure::LatticeModule  lattice;
	load_data( rootDir, loader, atomsUC, symmetry, kgrid, lattice);

	std::vector<int> symmetries = { 1, 0, 0,
							 0, 1, 0,
							 0, 0, 1,
							 0, 1, 0,
							 1, 0, 0,
							 0, 0, 1};

	std::vector<double> translations = { 0, 0, 0,
							 	 0, 0, 0};
	elephon::LatticeStructure::Symmetry overwrite;
	overwrite.initialize(1e-6,symmetries,translations,lattice,false);

	run_test( rootDir, loader, atomsUC, overwrite, kgrid, lattice, ref_forces_displ_0 );
}

BOOST_AUTO_TEST_CASE( Assemble_ForceConstantMatrix_full_symmetry )
{
	boost::filesystem::path rootDir = boost::filesystem::path(__FILE__).parent_path()
			/ "../data_for_testing/Al/vasp/phonon_test/";

	//reference data - pulled from the file - for irreducible displacement 0
	// which displaces atom 0 in direction x and is thus element mu2 = 0
	//Use it to compute the first column of phi explicitly
	std::vector<double> ref_forces_displ_0 =	//	Positions in unperturbed cell (mapped to the 1 UC):
	{-0.02181663 ,  0.00000000 , -0.00000000, // 0.00000000  0.00000000  0.00000000  R= 0 0 0 Atom : 0
	  0.02968642 , -0.00000000 ,  0.00000000, // 0.00000000 -0.25000000 -0.50000000  R= 0 0 0 Atom : 1
	 -0.00393490 , -0.00000000 ,  0.00000000, //-0.25000000  0.00000000 -0.50000000  R= 0 0 0 Atom : 2
	 -0.00393490 , -0.00000000 ,  0.00000000};//-0.25000000 -0.25000000  0.00000000  R= 0 0 0 Atom : 3

	double magn = -0.00123800*8.07786000;

	for ( auto &f : ref_forces_displ_0)
		f /= magn;

	std::shared_ptr<elephon::IOMethods::VASPInterface> loader;
	std::vector<elephon::LatticeStructure::Atom> atomsUC;
	elephon::LatticeStructure::Symmetry symmetry;
	elephon::LatticeStructure::RegularGrid kgrid;
	elephon::LatticeStructure::LatticeModule  lattice;
	load_data( rootDir, loader, atomsUC, symmetry, kgrid, lattice);

	run_test( rootDir, loader, atomsUC, symmetry, kgrid, lattice, ref_forces_displ_0 );
}
