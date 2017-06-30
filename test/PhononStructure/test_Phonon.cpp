/*	This file test_Phonon.cpp is part of elephon.
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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Input_test
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "fixtures/FixtureForceConstant.h"
#include "fixtures/DataLoader.h"
#include "fixtures/MockStartup.h"
#include "PhononStructure/Phonon.h"
#include <fstream>

BOOST_AUTO_TEST_CASE( Al_phonon_bands_gamma )
{
	test::fixtures::FixtureForceConstant ffc;
	auto fc = ffc.compute_fc_for_Al_gamma();

	std::vector<double> masses = {26.9815385};

	elephon::PhononStructure::Phonon ph;
	ph.initialize( fc, masses );

	BOOST_REQUIRE_EQUAL( ph.get_num_modes() , 3 );

	//compute 51 q points from 0,0,0 to 0.5,0,0
	int nq = 51;
	std::vector<double> qVect( nq * 3, 0.0 );
	for ( int i = 0 ; i < nq ; ++i)
		qVect[i*3] = 0.5*double(i)/double(nq-1);

	//load reference data from phonopy
	test::fixtures::MockStartup ms;
	auto testd = ms.get_data_for_testing_dir() / "Al" / "vasp" / "fcc_primitive" ;
	std::ifstream file( (testd / "phonopy_ref_phonons.dat").c_str() );
	if ( ! file.good() )
		throw std::runtime_error( (testd / "phonopy_ref_phonons.dat").string() + ": Could to be opened for reading");
	std::string line;
	double data[6];
	std::vector<double> refData( nq * 3, 0.0 );
	for ( int i = 0 ; i < nq; ++i)
	{
		std::getline(file,line);
		std::stringstream ss(line);
		ss >> data[0] >> data[1] >> data[2] >> data[3] >> data[4] >> data[5];
		if ( (std::abs(qVect[i*3]-data[0]) + std::abs(qVect[i*3+1]-data[1]) + std::abs(qVect[i*3+2]-data[2]) ) > 1e-6 )
			throw std::runtime_error( (testd / "phonopy_ref_phonons.dat").string() + ": q-vector differs from test");
		for ( int j = 0 ; j < 3 ; ++j)
			refData[i*3+j] = data[3+j];
	}

	std::vector<double> w;
	std::vector< std::complex<double> > eigenmodes;
	ph.compute_at_q( qVect, w, eigenmodes );
	assert( w.size() == ph.get_num_modes()*nq );
	for ( int iq = 0; iq < nq; ++iq)
		for ( int mu = 0; mu < ph.get_num_modes() ; ++mu)
			BOOST_CHECK_SMALL(refData[iq*ph.get_num_modes()+mu]-w[iq*ph.get_num_modes()+mu],0.1);
}
