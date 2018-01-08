/*	This file test_GradientFFTReciprocalGrid.cpp is part of elephon.
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
 *  Created on: Apr 30, 2017
 *      Author: A. Linscheid
 */
#include <boost/test/unit_test.hpp>
#include "ElectronicStructure/GradientFFTReciprocalGrid.h"
#include <assert.h>
#include <vector>

BOOST_AUTO_TEST_SUITE( GradientFFTReciprocalGrid )

BOOST_AUTO_TEST_CASE( Gradient_Cos_Sin )
{
	//Test the gradient of a system with cosine and sine data in two bands
	std::vector<int> grid({101,101,101});
	int gridnum = grid[0]*grid[1]*grid[2];
	std::vector<double> data(gridnum*2);

	for ( int k = 0 ; k < grid[2]; ++k)
		for ( int j = 0 ; j < grid[1]; ++j)
			for ( int i = 0 ; i < grid[0]; ++i)
			{
				data[((k*grid[1]+j)*grid[0]+i)*2+0] = -std::cos(4*M_PI*i/double(grid[0]))/(2.0);
				data[((k*grid[1]+j)*grid[0]+i)*2+1] = std::sin(4*M_PI*k/double(grid[2]))/(2.0);
			}

	std::vector<double> latticeMatrix = {   1.000000 , 0.000000 , 0.000000 ,
											0.000000 , 1.000000 , 0.000000 ,
											0.000000 , 0.000000 , 1.000000 };

	elephon::ElectronicStructure::GradientFFTReciprocalGrid gradE;
	gradE.compute_gradient(
			grid,
			latticeMatrix,
			2,
			data);

	BOOST_REQUIRE( gradE.get_data().size() == 2*gridnum*3 );

	std::vector<double> diff(6,0.0);
	for ( int k = 0 ; k < grid[2]; ++k)
		for ( int j = 0 ; j < grid[1]; ++j)
			for ( int i = 0 ; i < grid[0]; ++i)
			{
				//band 1
				int bandI = ((k*grid[1]+j)*grid[0]+i)*2+0;
				diff[0] += std::fabs(gradE.get_data()[bandI*3+0]-std::sin(4*M_PI*i/double(grid[0])));
				diff[1] += std::fabs(gradE.get_data()[bandI*3+1]);
				diff[2] += std::fabs(gradE.get_data()[bandI*3+2]);

				//band 2
				bandI = ((k*grid[1]+j)*grid[0]+i)*2+1;
				diff[3] += std::fabs(gradE.get_data()[bandI*3+0]);
				diff[4] += std::fabs(gradE.get_data()[bandI*3+1]);
				diff[5] += std::fabs(gradE.get_data()[bandI*3+2]-std::cos(4*M_PI*k/double(grid[2])));
			}

	double sum = 0;
	for ( auto a : diff )
		sum += a/double(gridnum);
	BOOST_REQUIRE( sum < 1e-6 );
}

BOOST_AUTO_TEST_CASE( Gradient_Cos_BndModel )
{
	//Test the gradient of a system with a cosine band model
	double const W1 = 200; // electron
	double const W2 = 100; // hole
	double const Ee = -50;
	double const Eg = 10;
	int nBnd = 2;
	std::vector<int> grid({101,100,100});
	int gridnum = grid[0]*grid[1]*grid[2];
	std::vector<double> data(gridnum*nBnd);

	for ( int k = 0 ; k < grid[2]; ++k)
		for ( int j = 0 ; j < grid[1]; ++j)
			for ( int i = 0 ; i < grid[0]; ++i)
			{
				double kx = (i<grid[0]/2?double(i):double(i)-grid[0])/double(grid[0]);
				double ky = (j<grid[1]/2?double(j):double(j)-grid[1])/double(grid[1]);
				  data[((k*grid[1]+j)*grid[0]+i)*nBnd+0] =
						  W1*(std::cos( 2*M_PI*kx )+std::cos( 2*M_PI*ky )-2.0)+Ee;
				  data[((k*grid[1]+j)*grid[0]+i)*nBnd+1] =
						  W2*(std::cos( 2*M_PI*kx )+std::cos( 2*M_PI*ky )-2.0)+Eg;
			}

	std::vector<double> latticeMatrix = {   1.000000 , 0.000000 , 0.000000 ,
											0.000000 , 1.000000 , 0.000000 ,
											0.000000 , 0.000000 , 1.000000 };

	elephon::ElectronicStructure::GradientFFTReciprocalGrid gradE;
	gradE.compute_gradient(
			grid,
			latticeMatrix,
			nBnd,
			data);

	BOOST_REQUIRE( gradE.get_data().size() == nBnd*gridnum*3 );

	std::vector<double> diff(6,0.0);
	for ( int k = 0 ; k < grid[2]; ++k)
		for ( int j = 0 ; j < grid[1]; ++j)
			for ( int i = 0 ; i < grid[0]; ++i)
			{
				double kx = (i<grid[0]/2?double(i):double(i)-grid[0])/double(grid[0]);
				double ky = (j<grid[1]/2?double(j):double(j)-grid[1])/double(grid[1]);

				//band 1
				int bandI = ((k*grid[1]+j)*grid[0]+i)*nBnd+0;
				diff[0] += std::fabs(gradE.get_data()[bandI*3+0]-
						  -W1*std::sin( 2*M_PI*kx ) );
				diff[1] += std::fabs(gradE.get_data()[bandI*3+1]-
						  -W1*std::sin( 2*M_PI*ky ));
				diff[2] += std::fabs(gradE.get_data()[bandI*3+2]);

				//band 2
				bandI = ((k*grid[1]+j)*grid[0]+i)*2+1;
				diff[3] += std::fabs(gradE.get_data()[bandI*3+0]-
						-W2*std::sin( 2*M_PI*kx ) );
				diff[4] += std::fabs(gradE.get_data()[bandI*3+1]-
						-W2*std::sin( 2*M_PI*ky ) );
				diff[5] += std::fabs(gradE.get_data()[bandI*3+2]);
			}

	double sum = 0;
	for ( auto a : diff )
		sum += a/double(gridnum);
	BOOST_REQUIRE_SMALL( sum , 1e-5 );
}

BOOST_AUTO_TEST_SUITE_END()
