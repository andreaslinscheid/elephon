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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Input_test
#include <boost/test/unit_test.hpp>
#include "ElectronicStructure/GradientFFTReciprocalGrid.h"
#include <assert.h>
#include <vector>

BOOST_AUTO_TEST_CASE( Gradient_Cos_Sin )
{
	//Test the gradient of a system with cosine and sine data in two bands
	std::vector<size_t> grid({101,101,101});
	size_t gridnum = grid[0]*grid[1]*grid[2];
	std::vector<double> data(gridnum*2);

	for ( size_t i = 0 ; i < grid[0]; ++i)
		for ( size_t j = 0 ; j < grid[1]; ++j)
			for ( size_t k = 0 ; k < grid[2]; ++k)
			{
				data[((i*grid[1]+j)*grid[2]+k)*2+0] = -std::cos(4*M_PI*i/double(grid[0]))/(4*M_PI);
				data[((i*grid[1]+j)*grid[2]+k)*2+1] = std::sin(4*M_PI*k/double(grid[2]))/(4*M_PI);
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
	for ( size_t i = 0 ; i < grid[0]; ++i)
		for ( size_t j = 0 ; j < grid[1]; ++j)
			for ( size_t k = 0 ; k < grid[2]; ++k)
			{
				//band 1
				size_t bandI = ((i*grid[1]+j)*grid[2]+k)*2+0;
				diff[0] += std::fabs(gradE.get_data()[bandI*3+0]-std::sin(4*M_PI*i/double(grid[0])));
				diff[1] += std::fabs(gradE.get_data()[bandI*3+1]);
				diff[2] += std::fabs(gradE.get_data()[bandI*3+2]);

				//band 2
				bandI = ((i*grid[1]+j)*grid[2]+k)*2+1;
				diff[3] += std::fabs(gradE.get_data()[bandI*3+0]);
				diff[4] += std::fabs(gradE.get_data()[bandI*3+1]);
				diff[5] += std::fabs(gradE.get_data()[bandI*3+2]-std::cos(4*M_PI*k/double(grid[2])));
			}

	double sum = 0;
	for ( auto a : diff )
		sum += a/double(gridnum);
	BOOST_REQUIRE( sum < 1e-6 );
}
