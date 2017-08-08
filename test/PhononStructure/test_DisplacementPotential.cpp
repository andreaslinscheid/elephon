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


BOOST_AUTO_TEST_CASE( build_Al_fcc_primitive )
{
	test::fixtures::FixtureForceConstant ff;
	elephon::PhononStructure::DisplacementPotential dvscf
		= ff.build_displ_pot_Al_fcc_primitive_vasp_sc2x2x2();

	BOOST_REQUIRE( dvscf.get_num_modes() == 3 );

	BOOST_REQUIRE( dvscf.get_num_R() == 2*2*2 );

	//TODO invent further checks that this is correct ...
}

BOOST_AUTO_TEST_CASE( plot_dvscf_q_Al_fcc_primitive )
{
	test::fixtures::MockStartup ms;
	auto rootDir = ms.get_data_for_testing_dir() / "Al" / "vasp" / "fcc_primitive" ;

	test::fixtures::FixtureForceConstant ff;
	elephon::PhononStructure::DisplacementPotential dvscf
		= ff.build_displ_pot_Al_fcc_primitive_vasp_sc2x2x2();

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
