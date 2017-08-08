/*	This file test_FFTInterface.cpp is part of elephon.
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
 *  Created on: Jul 8, 2017
 *      Author: A. Linscheid
 */
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Algorithms
#include <boost/test/unit_test.hpp>
#include "Algorithms/FFTInterface.h"
#include <vector>
#include <complex>

BOOST_AUTO_TEST_CASE( Cosine_tests_1D )
{

	elephon::Algorithms::FFTInterface fft;

	int nGrid = 100;
	int nBands = 1;
	std::vector<double> data(nGrid*nBands);
	for ( int i = 0 ; i < nGrid; ++i )
		data[i*nBands + 0] = std::cos( 2*M_PI* double(i<nGrid/2 ? i : i-nGrid)/double(nGrid) );

	std::vector< std::complex<double> > result;
	std::vector<int> dims{nGrid};
	fft.fft_data( dims, data, result, nBands, -1, true, 10 );

	BOOST_REQUIRE_EQUAL( result.size() , nGrid*nBands);
	double diff = 0;
	for ( int i = 0 ; i < nGrid; ++i )
		diff += std::abs(result[i*nBands + 0] - ((i == 1)||(i == nGrid-1) ? nGrid/2 : 0.0));

	BOOST_CHECK_SMALL( diff , 1e-6);
}

BOOST_AUTO_TEST_CASE( Cosine_sine_2Band_tests_1D )
{

	elephon::Algorithms::FFTInterface fft;

	int nGrid = 100;
	int nBands = 2;
	std::vector<double> data(nGrid*nBands);
	for ( int i = 0 ; i < nGrid; ++i )
	{
		data[i*nBands + 0] = std::cos( 2*M_PI* double(i<nGrid/2 ? i : i-nGrid)/double(nGrid) );
		data[i*nBands + 1] = std::sin( 2*M_PI* double(i<nGrid/2 ? i : i-nGrid)/double(nGrid) );
	}

	std::vector< std::complex<double> > result;
	std::vector<int> dims{nGrid};
	fft.fft_data( dims, data, result, nBands, -1, true, 10 );
	BOOST_REQUIRE_EQUAL( result.size() , nGrid*nBands);

	double diff = 0;
	for ( int i = 0 ; i < nGrid; ++i )
	{
		diff += std::abs(result[i*nBands + 0] - ((i == 1)||(i == nGrid-1) ? nGrid/2 : 0.0));
		diff += std::abs(result[i*nBands + 1] - std::complex<double>(0,-(i == 1? nGrid/2 : 0.0)+(i == nGrid-1?nGrid/2 : 0.0)) );
	}

	BOOST_CHECK_SMALL( diff , 1e-6);
}

BOOST_AUTO_TEST_CASE( Cosine_sine_2Band_tests_2D )
{

	elephon::Algorithms::FFTInterface fft;

	std::vector<int> grid({101,100});
	int nBands = 2;
	int nGrid = grid[0]*grid[1];
	std::vector<double> data(nGrid*nBands);

	//Column major data arrangement
	for ( int j = 0 ; j < grid[1]; ++j)
		for ( int i = 0 ; i < grid[0]; ++i)
		{
			data[(j*grid[0]+i)*nBands+0] = std::cos(2*M_PI*i/double(grid[0]));
			data[(j*grid[0]+i)*nBands+1] = std::sin(2*M_PI*j/double(grid[1]));
		}

	std::vector< std::complex<double> > result;
	fft.fft_data( grid, data, result, nBands, -1, false, 10 );

	BOOST_REQUIRE_EQUAL( result.size() , nGrid*nBands);

	double diff = 0;
	for ( int j = 0 ; j < grid[1]; ++j)
		for ( int i = 0 ; i < grid[0]; ++i)
		{
			diff += std::abs(result[(j*grid[0]+i)*nBands + 0] -
					(j==0?((i == 1)||(i == grid[0]-1) ? nGrid/2 : 0.0):0.0));
			diff += std::abs(result[(j*grid[0]+i)*nBands + 1] -
					(i==0?(std::complex<double>(0,-(j == 1? nGrid/2 : 0.0)+(j == grid[1]-1?nGrid/2 : 0.0))):0.0) );
		}

	BOOST_CHECK_SMALL( diff , 1e-6);
}

BOOST_AUTO_TEST_CASE( sparse_data_test )
{
	elephon::Algorithms::FFTInterface fft;

	int nBands = 1;
	std::vector<double> data(2*nBands);

	std::vector<int> grid({101,100,150});

	//Create the Fourier transform of cos(x) in 3D
	std::vector<int> fftmap = {
			1			, 0,0,
			grid[0]-1	, 0,0 };

	data[0] = 0.5;
	data[1] = 0.5;

	std::vector< std::complex<double> > result;
	fft.fft_sparse_data( fftmap, data, nBands, -1, result, grid, false, 10 );

	int ngrid = grid[0]*grid[1]*grid[2];
	BOOST_REQUIRE_EQUAL( result.size() , ngrid*nBands);

	double diff = 0;
	for ( int k = 0 ; k < grid[2]; ++k)
		for ( int j = 0 ; j < grid[1]; ++j)
			for ( int i = 0 ; i < grid[0]; ++i)
				diff += std::abs(result[((k*grid[1]+j)*grid[0]+i)*nBands + 0] - std::cos(2*M_PI*i/double(grid[0])));

	BOOST_CHECK_SMALL( diff , 1e-6);
}
