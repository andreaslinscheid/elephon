/*	This file test_ReadVASPWaveFunction.cpp is part of elephon.
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
 *  Created on: May 18, 2017
 *      Author: A. Linscheid
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Input_test
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "IOMethods/ReadVASPWaveFunction.h"
#include <vector>
#include <complex>

BOOST_AUTO_TEST_CASE( Read_VASP_Al_wavefunctions )
{
	using namespace boost::filesystem;
	path p(__FILE__);
	path dir = p.parent_path();
	path test_Al_wfctfile = dir / "WAVECAR_Al_test.dat";

	elephon::IOMethods::ReadVASPWaveFunction wfctreader;
	wfctreader.prepare_wavecar( test_Al_wfctfile.string() );

	int nBnd = 10;
	BOOST_REQUIRE( wfctreader.get_num_bands() == nBnd );
	BOOST_REQUIRE( wfctreader.get_num_spins() == 1 );
	BOOST_REQUIRE( wfctreader.get_num_kpts() == 364 );

	std::ifstream compareFile( (dir / "wavecar_cmp.dat").c_str() );

	//Check if the energy was read in correctly
	std::string buffer;
	std::getline( compareFile, buffer );
	//first is the energy of band 0 and k point 4
	double energy_k4_b0;
	int npw =  1788;
	std::stringstream ss1(buffer);
	ss1 >> npw >> energy_k4_b0;
	BOOST_REQUIRE( std::abs(energy_k4_b0 - wfctreader.get_energies()[4*nBnd+0]) < 1e-10 );

	//Check if the k-point was read in correctly
	double k[3];
	std::getline( compareFile, buffer );
	std::stringstream ss2(buffer);
	ss2 >> k[0] >> k[1] >> k[2];
	BOOST_REQUIRE( (std::abs( k[0] - wfctreader.get_k_points()[4*3+0]) < 1e-10)
			&& (std::abs( k[1] - wfctreader.get_k_points()[4*3+1]) < 1e-10)
			&& (std::abs( k[2] - wfctreader.get_k_points()[4*3+2]) < 1e-10));

	std::vector<int> kpts = {1,4};
	std::vector<int> bands = {0,1,4};
	std::vector< std::complex<float> > wfct;
	std::vector< std::vector<int> > fourier;
	std::vector<int> fft;
	wfctreader.read_wavefunction( kpts, bands, wfct, fourier, fft );

	//Check the second kpt we wanted to have which is ik=4
	BOOST_REQUIRE( static_cast<int>(fourier[1].size())/3 == npw );

	int nelemWfct = fft[0]*fft[1]*fft[2];

	//Read the reference data ...
	std::getline( compareFile, buffer );
	std::stringstream ss3(buffer);
	std::vector< float > pwcoeffr( (std::istream_iterator<float>(ss3)), std::istream_iterator<float>() );
	std::vector< std::complex<float> > pwcoeff(
			reinterpret_cast<std::complex<float> *>(pwcoeffr.data()),
			reinterpret_cast<std::complex<float> *>(pwcoeffr.data())+npw);

	//...and compare the second k point (ik=4) and the zeroth band (ibnd = 0)
	float diff = 0;
	for ( int i = 0 ; i < npw ; ++i)
		diff += std::abs(wfct[(1*bands.size()+0)*nelemWfct+i]-pwcoeff[i]);

	//Note on realistic requirements for equivalence: we compare (float )->(float) and (float)->(text)->(float).
	//The text is reasonably large so that we should roughly keep float precision but its not the _same_ number.
	float norm = 0;
	for ( auto n : pwcoeff )
		norm += std::abs(n);
	BOOST_REQUIRE( diff/norm < 1e-5 );

	std::getline( compareFile, buffer );
	std::stringstream ss4(buffer);
	std::vector< int > fourierMap( (std::istream_iterator<int>(ss4)), std::istream_iterator<int>() );
	int diffMap = 0;
	for ( int i = 0 ; i < npw*3 ; ++i)
		diffMap += std::abs( fourier[1][i] - (fourierMap[i]-1) );
		// the -1 above is the fortran style reference data that counts from 1, not zero

	BOOST_REQUIRE( diffMap == 0);
}
