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
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "fixtures/DataLoader.h"
#include "fixtures/MockStartup.h"
#include "PhononStructure/Phonon.h"
#include "fixtures/scenarios.h"
#include <fstream>

BOOST_AUTO_TEST_SUITE( Phonon )

BOOST_AUTO_TEST_CASE( Al_phonon_bands_gamma )
{
	auto res = elephon::test::fixtures::scenarios::load_Al_fcc_primitive_vasp_sc2x2x2();

	auto ph = res->get_phonon_obj();
	BOOST_REQUIRE_EQUAL( ph->get_num_modes() , 3 );

	//compute 51 q points from 0,0,0 to 0.5,0,0
	int nq = 51;
	std::vector<double> qVect( nq * 3, 0.0 );
	for ( int i = 0 ; i < nq ; ++i)
		qVect[i*3] = 0.5*double(i)/double(nq-1);

	//load reference data from phonopy
	elephon::test::fixtures::MockStartup ms;
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

	elephon::Auxillary::alignedvector::DV w;
	elephon::Auxillary::alignedvector::ZV eigenmodes;
	ph->compute_at_q( qVect, w, eigenmodes );
	assert( w.size() == ph->get_num_modes()*nq );
	for ( int iq = 0; iq < nq; ++iq)
	{
		for ( int mu = 0; mu < ph->get_num_modes() ; ++mu)
			BOOST_CHECK_SMALL(refData[iq*ph->get_num_modes()+mu]-w[iq*ph->get_num_modes()+mu],0.1);
		}

	for ( int iq = 0; iq < nq; ++iq)
	{
		std::vector<std::complex<double>> u1(3),u2(3),u3(3);
		for ( int i = 0; i < 3; ++i)
		{
			u1[i] = eigenmodes[iq*9+0*3+i];
			u2[i] = eigenmodes[iq*9+1*3+i];
			u3[i] = eigenmodes[iq*9+2*3+i];
		}

		BOOST_CHECK_CLOSE( std::abs(u1[0]*u1[0]+u1[1]*u1[1]+u1[2]*u1[2]) , 1.0, 0.02);
		BOOST_CHECK_CLOSE( std::abs(u2[0]*u2[0]+u2[1]*u2[1]+u2[2]*u2[2]) , 1.0, 0.02);
		BOOST_CHECK_CLOSE( std::abs(u3[0]*u3[0]+u3[1]*u3[1]+u3[2]*u3[2]) , 1.0, 0.02);

		BOOST_CHECK_SMALL( std::abs(u1[0]*u2[0]+u1[1]*u2[1]+u1[2]*u2[2]) , 0.0002);
		BOOST_CHECK_SMALL( std::abs(u1[0]*u3[0]+u1[1]*u3[1]+u1[2]*u3[2]) , 0.0002);
		BOOST_CHECK_SMALL( std::abs(u2[0]*u3[0]+u2[1]*u3[1]+u2[2]*u3[2]) , 0.0002);
	}
}

BOOST_AUTO_TEST_CASE( Al_phonon_derivative )
{
	auto res = elephon::test::fixtures::scenarios::load_Al_fcc_primitive_vasp_sc2x2x2();

	auto ph = res->get_phonon_obj();

	int nq = 1510;
	std::vector<double> qVect( nq * 3, 0.0 );
	for ( int i = 0 ; i < nq ; ++i)
		qVect[i*3] = 0.5*double(i)/double(nq-1);

	elephon::Auxillary::alignedvector::DV omega, domegadq;
	elephon::Auxillary::alignedvector::ZV em;
	ph->compute_at_q(qVect, omega, em);
	ph->evaluate_derivative(qVect, domegadq);

	std::vector<double> v1{0.5/double(nq-1), 0.0, 0.0};
	auto const & kgrid = res->get_electronic_bands_obj()->get_grid();
	kgrid.get_lattice().reci_direct_to_cartesian_2pibya(v1);

	double diff = 0.0;
	for ( int iq = 1; iq < nq-1; ++iq)
	{
		for ( int mu = 0; mu < ph->get_num_modes() ; ++mu)
		{
			double norm = std::sqrt(std::pow(v1[0],2)+std::pow(v1[1],2)+std::pow(v1[2],2));
			double directionDerv =
					(v1[0]*domegadq[(iq*ph->get_num_modes()+mu)*3]
					+ v1[1]*domegadq[(iq*ph->get_num_modes()+mu)*3+1]
					+ v1[2]*domegadq[(iq*ph->get_num_modes()+mu)*3+2])/norm;
			double simplyDeriv = (omega[(iq+1)*ph->get_num_modes()+mu] - omega[iq*ph->get_num_modes()+mu])
					/ norm;

			diff += std::abs(directionDerv - simplyDeriv) / (nq-2);
		}
	}

	BOOST_CHECK_SMALL(diff, 1e-2);
}

BOOST_AUTO_TEST_SUITE_END()
