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
#include <boost/test/unit_test.hpp>
#include "Algorithms/FFTInterface.h"
#include <vector>
#include <complex>

BOOST_AUTO_TEST_SUITE( FFTInterface )

BOOST_AUTO_TEST_CASE( fft_interpolate_cos_shift_1D )
{
	std::vector<int> grid{20};
	std::vector<double> data(grid[0]);
	std::vector<double> sIn {0.50};
	for ( int i = 0 ; i < grid[0]; ++i)
		data[i] = std::cos(2*M_PI*((i+sIn[0])/double(grid[0])));

	std::vector<int> gridOut{50};
	std::vector<double> sOut{0.0};
	int nG = gridOut[0];
	elephon::Algorithms::FFTInterface fft;
	fft.fft_interpolate(grid, sIn, data, gridOut, sOut, data, 1);

	double diff = 0;
	for ( int i = 0 ; i < gridOut[0]; ++i)
		diff += std::abs(data[i] - std::cos(2*M_PI*((i+sOut[0])/double(gridOut[0]))));

   BOOST_CHECK_SMALL( diff/nG , 1e-6);
}

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
	fft.plan_fft( dims, nBands, -1, true, 10);
	fft.fft_data(data, result, -1);

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
	fft.plan_fft(dims, nBands, -1, true, 10);
	fft.fft_data(data, result, -1);
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
	fft.plan_fft(grid, nBands, -1, false, 10);
	fft.fft_data(data, result, -1);

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
	fft.plan_fft(grid, nBands, -1, false, 10);
	fft.fft_sparse_data( fftmap, grid, data, -1, result);

	int ngrid = grid[0]*grid[1]*grid[2];
	BOOST_REQUIRE_EQUAL( result.size() , ngrid*nBands);

	double diff = 0;
	for ( int k = 0 ; k < grid[2]; ++k)
		for ( int j = 0 ; j < grid[1]; ++j)
			for ( int i = 0 ; i < grid[0]; ++i)
				diff += std::abs(result[((k*grid[1]+j)*grid[0]+i)*nBands + 0] - std::cos(2*M_PI*i/double(grid[0])));

	BOOST_CHECK_SMALL( diff , 1e-6);
}

BOOST_AUTO_TEST_CASE( fft_interpolate_constant )
{
       int nBands = 3;
       std::vector<int> grid({1,1,1});
       std::vector<double> data(grid[0]*grid[1]*grid[2]*nBands);
       data[0] = 1;
       data[1] = 2;
       data[2] = -2;

       std::vector<int> gridOut({3,2,4});
       int nG = gridOut[0]*gridOut[1]*gridOut[2];
       elephon::Algorithms::FFTInterface fft;
       auto nullv = std::vector<double>{0.0, 0.0, 0.0};
       fft.fft_interpolate(grid, nullv, data, gridOut, nullv, data, nBands);

       BOOST_REQUIRE_EQUAL(data.size(), nBands*nG);
       for ( int k = 0 ; k < gridOut[2]; ++k)
               for ( int j = 0 ; j < gridOut[1]; ++j)
                       for ( int i = 0 ; i < gridOut[0]; ++i)
                       {
                               int cnsq = i + gridOut[0]*(j + gridOut[1]*k);
                               BOOST_CHECK_SMALL(data[cnsq*nBands+0] - 1, 1e-6);
                               BOOST_CHECK_SMALL(data[cnsq*nBands+1] - 2, 1e-6);
                               BOOST_CHECK_SMALL(data[cnsq*nBands+2] + 2, 1e-6);
                       }
}

BOOST_AUTO_TEST_CASE( fft_interpolate_cos_sin )
{
	int nBands = 2;
	std::vector<int> grid({11,8,9});
	std::vector<double> data(grid[0]*grid[1]*grid[2]*nBands);
	for ( int k = 0 ; k < grid[2]; ++k)
		   for ( int j = 0 ; j < grid[1]; ++j)
				   for ( int i = 0 ; i < grid[0]; ++i)
				   {
						   int cnsq = i + grid[0]*(j + grid[1]*k);
						   data[cnsq*nBands+0] = std::cos(2*M_PI*i/double(grid[0]));
						   data[cnsq*nBands+1] = std::pow(std::sin(2*M_PI*k/double(grid[2])),2);
				   }

	std::vector<int> gridOut({15,7,28});
	int nG = gridOut[0]*gridOut[1]*gridOut[2];
	elephon::Algorithms::FFTInterface fft;
	std::vector<double> dataOut;
	fft.fft_interpolate(
				   grid,
				   std::vector<double>{0.0, 0.0, 0.0},
				   data,
				   gridOut,
				   std::vector<double>{0.0, 0.0, 0.0},
				   dataOut,
				   nBands);

	BOOST_REQUIRE_EQUAL(dataOut.size(), nBands*nG);

	double diff = 0;
	for ( int k = 0 ; k < gridOut[2]; ++k)
		for ( int j = 0 ; j < gridOut[1]; ++j)
			for ( int i = 0 ; i < gridOut[0]; ++i)
			{
				int cnsq = i + gridOut[0]*(j + gridOut[1]*k);
				diff += std::abs(dataOut[cnsq*nBands+0] - std::cos(2*M_PI*i/double(gridOut[0])));
				diff += std::abs(dataOut[cnsq*nBands+1] - std::pow(std::sin(2*M_PI*k/double(gridOut[2])),2));
			}

   BOOST_CHECK_SMALL( diff/nBands/nG , 1e-6);
}

BOOST_AUTO_TEST_CASE( fft_interpolate_cos_sin_shift )
{
	std::vector<int> grid({11,7,5});
	std::vector<double> data(grid[0]*grid[1]*grid[2]);
	std::vector<double> sIn {0.50, 0.75, 0.00};
	for ( int k = 0 ; k < grid[2]; ++k)
		for ( int j = 0 ; j < grid[1]; ++j)
			for ( int i = 0 ; i < grid[0]; ++i)
				data[i + grid[0]*(j + grid[1]*k)] = std::cos(2*M_PI*((i+sIn[0])/double(grid[0])))
													+std::cos(2*M_PI*((j+sIn[1])/double(grid[1])));

	std::vector<int> gridOut({15,27,23});
	std::vector<double> sOut{0.25, 0.00, 0.00};
	int nG = gridOut[0]*gridOut[1]*gridOut[2];
	elephon::Algorithms::FFTInterface fft;
	std::vector<double> dataOut;
	fft.fft_interpolate(grid, sIn, data, gridOut, sOut, dataOut, 1);

	BOOST_REQUIRE_EQUAL(dataOut.size(), nG);

	double diff = 0;
	for ( int k = 0 ; k < gridOut[2]; ++k)
		for ( int j = 0 ; j < gridOut[1]; ++j)
			for ( int i = 0 ; i < gridOut[0]; ++i)
				diff += std::abs(dataOut[i + gridOut[0]*(j + gridOut[1]*k)]
									  - std::cos(2*M_PI*((i+sOut[0])/double(gridOut[0])))
									  - std::cos(2*M_PI*((j+sOut[1])/double(gridOut[1]))));

   BOOST_CHECK_SMALL( diff/nG , 1e-6);
}

BOOST_AUTO_TEST_CASE( fft_hessian_1D_cos )
{
	std::vector<int> grid{250};
	int nG = grid[0];
	std::vector<double> data(nG);
	for ( int i = 0 ; i < grid[0]; ++i)
		data[i] = std::cos(2*M_PI*i/double(grid[0]));
	std::vector<double> latticeMatrix{1.0};
	elephon::Algorithms::FFTInterface fft;
	std::vector<double> dataOut;
	fft.fft_hessian(latticeMatrix, grid, data, dataOut, 1);
	double diff = 0;
	for ( int i = 0 ; i < grid[0]; ++i)
		diff += std::abs(dataOut[i]-(-std::cos(2*M_PI*i/double(grid[0]))));
	BOOST_CHECK_SMALL( diff/nG , 1e-6);
}

BOOST_AUTO_TEST_CASE( fft_hessian_2D_cos )
{
	std::vector<int> grid{128, 128};
	int nG = grid[0]*grid[1];
	int D = 2;
	std::vector<double> data(nG);
	for ( int j = 0 ; j < grid[1]; ++j)
		for ( int i = 0 ; i < grid[0]; ++i)
			data[i + grid[0]*j] = std::cos(2*M_PI*i/double(grid[0]))+std::cos(2*M_PI*j/double(grid[1]));
	std::vector<double> latticeMatrix{	1.0, 0.0,
										0.0, 1.0 };
	elephon::Algorithms::FFTInterface fft;
	std::vector<double> dataOut;
	fft.fft_hessian(latticeMatrix, grid, data, dataOut, 1);
	// compare against analytic formula in band 0
	double diff = 0;
	for ( int j = 0 ; j < grid[1]; ++j)
		for ( int i = 0 ; i < grid[0]; ++i)
			for ( int jgr = 0 ; jgr < D; ++jgr)
				for ( int igr = 0 ; igr < D; ++igr)
					diff += std::abs(dataOut[igr+D*(jgr+D*(i + grid[0]*j))] -  (
							((igr == 0)&&(jgr == 0) ? - std::cos((2*M_PI*i)/grid[0]) : 0.0) +
							((igr == 1)&&(jgr == 1) ? - std::cos((2*M_PI*j)/grid[1]) : 0.0)  ) );
	BOOST_CHECK_SMALL( diff/nG/std::pow(D,2) , 1e-6);
}

BOOST_AUTO_TEST_CASE( fft_hessian_3D_cos )
{
	std::vector<int> grid{128, 128, 64};
	int nG = grid[0]*grid[1]*grid[2];
	int nB = 2;
	int D = 3;
	std::vector<double> data(nG*nB);
	for ( int k = 0 ; k < grid[2]; ++k)
		for ( int j = 0 ; j < grid[1]; ++j)
			for ( int i = 0 ; i < grid[0]; ++i)
			{
				int cnsq = i + grid[0]*(j + grid[1]*k);
				data[cnsq*nB+0] = std::cos(2*M_PI*i/double(grid[0]))+std::cos(2*M_PI*j/double(grid[1]));
				data[cnsq*nB+1] = std::pow(std::sin(2*M_PI*k/double(grid[2])),2);
			}
	std::vector<double> latticeMatrix{	1.0, 0.0, 0.0,
										0.0, 1.0, 0.0,
										0.0, 0.0, 1.0 };
	elephon::Algorithms::FFTInterface fft;
	std::vector<double> dataOut;
	fft.fft_hessian(latticeMatrix, grid, data, dataOut, nB);
	double diff = 0;
	for ( int k = 0 ; k < grid[2]; ++k)
		for ( int j = 0 ; j < grid[1]; ++j)
			for ( int i = 0 ; i < grid[0]; ++i)
			{
				int cnsq = i + grid[0]*(j + grid[1]*k);
				// compare against analytic formula in band 0
				for ( int jgr = 0 ; jgr < D; ++jgr)
					for ( int igr = 0 ; igr < D; ++igr)
						diff += std::abs(dataOut[igr+D*(jgr+D*(0+nB*cnsq))] - (
								((igr == 0)&&(jgr == 0) ? - std::cos((2*M_PI*i)/grid[0]) : 0.0) +
								((igr == 1)&&(jgr == 1) ? - std::cos((2*M_PI*j)/grid[1]) : 0.0)  ) );
				// and in band 1
				for ( int jgr = 0 ; jgr < D; ++jgr)
					for ( int igr = 0 ; igr < D; ++igr)
						diff += std::abs(dataOut[igr+D*(jgr+D*(1+nB*cnsq))] - ((igr == 2)&&(jgr == 2) ? 1.0 : 0.0 )*
								2.0*(std::pow(std::cos((2*M_PI*k)/grid[2]),2) - std::pow(std::sin((2*M_PI*k)/grid[2]),2)));
			}
	BOOST_CHECK_SMALL( diff/nG/nB/std::pow(D,2) , 1e-6);
}

BOOST_AUTO_TEST_CASE( fft_hessian_3D_cos_non_cubic_cell )
{
	std::vector<double> latticeMatrix{	1.000000, 1.000000, 0.000000,
									   -1.000000, 1.000000, 0.000000,
										0.000000, 0.000000, 1.000000 };
	elephon::LatticeStructure::LatticeModule latticeMod( latticeMatrix );

	std::vector<int> grid{64, 64, 64};
	int nG = grid[0]*grid[1]*grid[2];
	int D = 3;
	std::vector<double> data(nG);

	// here we initialize data in the cubic cell. In terms of the lattice basis
	// this will of cause be tilted.
	std::vector<double> cubicVect(D);
	for ( int k = 0 ; k < grid[2]; ++k)
		for ( int j = 0 ; j < grid[1]; ++j)
			for ( int i = 0 ; i < grid[0]; ++i)
			{
				cubicVect = {i/double(grid[0]), j/double(grid[1]), k/double(grid[2])};
				latticeMod.direct_to_cartesian_angstroem(cubicVect);
				int cnsq = i + grid[0]*(j + grid[1]*k);
				data[cnsq] = std::cos(2*M_PI*cubicVect[0])+std::cos(2*M_PI*cubicVect[1])
								+std::cos(2*M_PI*cubicVect[2]);
			}

	auto B = latticeMod.get_reciprocal_latticeMatrix();
	for ( auto &b : B)
		b /= latticeMod.get_alat();

	elephon::Algorithms::FFTInterface fft;
	std::vector<double> dataOut;
	fft.fft_hessian(B, grid, data, dataOut, 1);

	// independent on the coordinate system, this must give the correct derivative
	double diff = 0;
	for ( int k = 0 ; k < grid[2]; ++k)
		for ( int j = 0 ; j < grid[1]; ++j)
			for ( int i = 0 ; i < grid[0]; ++i)
			{
				cubicVect = {i/double(grid[0]), j/double(grid[1]), k/double(grid[2])};
				latticeMod.direct_to_cartesian_angstroem(cubicVect);
				int cnsq = i + grid[0]*(j + grid[1]*k);
				double x = 2*M_PI*cubicVect[0];
				double y = 2*M_PI*cubicVect[1];
				double z = 2*M_PI*cubicVect[2];
				for ( int jgr = 0 ; jgr < D; ++jgr)
					for ( int igr = 0 ; igr < D; ++igr)
						diff += std::abs(dataOut[igr+D*(jgr+D*cnsq)] - (
								((igr == 0)&&(jgr == 0) ? - std::cos(x) : 0.0) +
								((igr == 1)&&(jgr == 1) ? - std::cos(y) : 0.0) +
								((igr == 2)&&(jgr == 2) ? - std::cos(z) : 0.0)  ) );
			}
	BOOST_CHECK_SMALL( diff/nG/std::pow(D,2) , 1e-6);
}

BOOST_AUTO_TEST_CASE( sparse_data_parallel_test )
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

	fft.plan_fft(grid, nBands, -1, false, 10);

	double diff = 0;
	#pragma omp parallel
	{
		// thread private result
		std::vector< std::complex<double> > result;
		double diff_thread_local = 0;

		fft.fft_sparse_data( fftmap, grid, data, -1, result);

		int ngrid = grid[0]*grid[1]*grid[2];
		BOOST_REQUIRE_EQUAL( result.size() , ngrid*nBands);

		for ( int k = 0 ; k < grid[2]; ++k)
			for ( int j = 0 ; j < grid[1]; ++j)
				for ( int i = 0 ; i < grid[0]; ++i)
					diff_thread_local += std::abs(result[((k*grid[1]+j)*grid[0]+i)*nBands + 0] - std::cos(2*M_PI*i/double(grid[0])));

		#pragma omp atomic
		diff += diff_thread_local;
	}

	BOOST_CHECK_SMALL( diff , 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()
