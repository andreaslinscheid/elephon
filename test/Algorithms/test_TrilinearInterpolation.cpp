
/*	This file test_TrilinearInterpolation.cpp is part of elephon.
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
 *  Created on: Apr 28, 2017
 *      Author: A. Linscheid
 */
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Input_test
#include <boost/test/unit_test.hpp>
#include "Algorithms/TrilinearInterpolation.h"
#include <LatticeStructure/RegularSymmetricGrid.h>
#include "LatticeStructure/Symmetry.h"
#include "LatticeStructure/LatticeModule.h"
#include <vector>
#include <iostream>

BOOST_AUTO_TEST_CASE( Linear_Interpolation_Exact )
{
	const double accuracy = 10e-8;
	//Test that our linear interpolation interpolates
	std::vector<int> testGrid = {2,2,2};

	std::vector<double> testPoints = {  0.0,0.0,0.0 ,
										0.1,0.05,0.025 };

	elephon::LatticeStructure::RegularBareGrid regularGrid;
	regularGrid.initialize( testGrid );
	elephon::Algorithms::TrilinearInterpolation triLin(regularGrid);
	std::vector<int> indices;
	triLin.data_query(testPoints,indices);

	//Should be all in the first cell
	BOOST_REQUIRE_EQUAL( indices.size() , 8 );

	//generate linear data in x
	std::vector<double> data(testGrid[0]*testGrid[1]*testGrid[2]);
	for ( int i = 0 ; i < testGrid[2]; ++i)
		for ( int j = 0 ; j < testGrid[1]; ++j)
			for ( int k = 0 ; k < testGrid[0]; ++k)
				data[(i*testGrid[1]+j)*testGrid[0]+k] = double(k)/double(testGrid[0]);

	//copy the necessary indices
	std::vector<double> neededData, interpolatedData;
	for ( auto i : indices )
		neededData.push_back(data[i]);

	triLin.interpolate(1,neededData,interpolatedData);
	BOOST_REQUIRE_EQUAL( interpolatedData.size() , 2 );
	BOOST_REQUIRE( std::fabs(interpolatedData[0] - 0.0) < accuracy );
	BOOST_REQUIRE( std::fabs(interpolatedData[1] - 0.1) < accuracy );

	//repeat for y - note the second point is also y=0.05
	for ( int i = 0 ; i < testGrid[2]; ++i)
		for ( int j = 0 ; j < testGrid[1]; ++j)
			for ( int k = 0 ; k < testGrid[0]; ++k)
				data[(i*testGrid[1]+j)*testGrid[0]+k] = double(j)/double(testGrid[1]);
	neededData.clear();
	for ( auto i : indices )
		neededData.push_back(data[i]);
	triLin.interpolate(1,neededData,interpolatedData);
	BOOST_REQUIRE( std::fabs(interpolatedData[0] - 0.0) < accuracy );
	BOOST_REQUIRE( std::fabs(interpolatedData[1] - 0.05) < accuracy );

	//repeat for z - note the second point is also z=0.025
	for ( int i = 0 ; i < testGrid[2]; ++i)
		for ( int j = 0 ; j < testGrid[1]; ++j)
			for ( int k = 0 ; k < testGrid[0]; ++k)
				data[(i*testGrid[1]+j)*testGrid[0]+k] = double(i)/double(testGrid[2]);
	neededData.clear();
	for ( auto i : indices )
		neededData.push_back(data[i]);
	triLin.interpolate(1,neededData,interpolatedData);
	BOOST_REQUIRE( std::fabs(interpolatedData[0] - 0.0) < accuracy );
	BOOST_REQUIRE( std::fabs(interpolatedData[1] - 0.025) < accuracy );
}

BOOST_AUTO_TEST_CASE( Linear_Interpolation_Blockdata )
{
	//Here we test the linear interpolation of several data values per grid point
	//We also use awkward grid dimension to make sure this works too.
	const double accuracy = 10e-8;
	//Test that our linear interpolation interpolates
	std::vector<int> testGrid = {7,5,2};

	std::vector<double> testPoints = {   0.0, 0.0,  0.0 ,//p1
										 0.1, 0.05, 0.025,//p2
										-0.5,-0.49, 0.7}; //p3

	elephon::LatticeStructure::RegularBareGrid regularGrid;
	regularGrid.initialize( testGrid );
	elephon::Algorithms::TrilinearInterpolation triLin(regularGrid);
	std::vector<int> indices;
	triLin.data_query( testPoints,indices);

	//p1 and p2 be all in the cell starting with the vector having the indices (0,0,0)
	// and p3 ~= 0.5,0.51,0.7 should be in the cell starting with the vector having indices (3,2,1)
	// thus, they share no points. The total number of points should be thus 8*2=16
	BOOST_REQUIRE_EQUAL( indices.size() , 16 );

	int nBlock = 2;
	//generate linear data in x,y,z
	std::vector<double> data(testGrid[0]*testGrid[1]*testGrid[2]*nBlock);
	for ( int i = 0 ; i < testGrid[2]; ++i)
		for ( int j = 0 ; j < testGrid[1]; ++j)
			for ( int k = 0 ; k < testGrid[0]; ++k)
				for ( int ib = 0 ; ib < nBlock; ++ib)
					data[((i*testGrid[1]+j)*testGrid[0]+k)*nBlock + ib] =
							double(ib)+
						double(k)/double(testGrid[0]);

	//copy the necessary indices
	std::vector<double> neededData, interpolatedData;
	for ( auto i : indices )
		for ( int ib = 0 ; ib < nBlock; ++ib)
			neededData.push_back(data[i*nBlock+ib]);
	triLin.interpolate(nBlock,neededData,interpolatedData);
	BOOST_REQUIRE( std::fabs(interpolatedData[0] - 0.0) < accuracy );
	BOOST_REQUIRE( std::fabs(interpolatedData[1] - 1.0) < accuracy );
	BOOST_REQUIRE( std::fabs(interpolatedData[2] - 0.1) < accuracy );
	BOOST_REQUIRE( std::fabs(interpolatedData[3] - 1.1) < accuracy );
	BOOST_REQUIRE( std::fabs(interpolatedData[4] - 0.5) < accuracy );
	BOOST_REQUIRE( std::fabs(interpolatedData[5] - 1.5) < accuracy );
}

BOOST_AUTO_TEST_CASE( Linear_Interpolation_IntegralTest )
{
	//Here we use a dense grid and sample a cosine function.
	//Then we use an even denser grid and use the trilinear interpolation.
	//Given that the first function is already dense, the integral of the linearly interpolated function must
	//not be very different.
	//In the first block, we use |cos(2pi x)||cos(2pi y)||cos(2pi z)|. The integral of xi=0,0.5 gives 1/(pi^3)
	//In the second block, we use sin(2pi x)sin(2pi y)sin(2pi z)
	//In the third block we use 1. The integral is obviously 1/8

	const double accuracy = 1e-4;

	std::vector<int> testGrid = {60,60,60};

	//Generate a denser mesh of 150^3 points in the first 1/2,1/2,1/3 wedge of the unit cell
	std::vector<double> testPoints(testGrid[0]*testGrid[1]*testGrid[2] *3);

	int nBlock = 3;
	std::vector<double> data(testGrid[0]*testGrid[1]*testGrid[2]*nBlock);
	double integralBlock1 = 0;
	double integralBlock2 = 0;
	double integralBlock3 = 0;
	double dV = 1.0/double(testGrid[0]*testGrid[1]*testGrid[2]);
	for ( int i = 0 ; i < testGrid[0]; ++i)
		for ( int j = 0 ; j < testGrid[1]; ++j)
			for ( int k = 0 ; k < testGrid[2]; ++k)
			{
				double x = double(i)/double(testGrid[0]);
				double y = double(j)/double(testGrid[1]);
				double z = double(k)/double(testGrid[2]);
				data[((i*testGrid[1]+j)*testGrid[0]+k)*nBlock + 0] =
						std::fabs(std::cos(2.0*M_PI*x))*
						std::fabs(std::cos(2.0*M_PI*y))*
						std::fabs(std::cos(2.0*M_PI*z));
				data[((i*testGrid[1]+j)*testGrid[0]+k)*nBlock + 1] =
						std::sin(2.0*M_PI*x)*
						std::sin(2.0*M_PI*y)*
						std::sin(2.0*M_PI*z);
				data[((i*testGrid[1]+j)*testGrid[0]+k)*nBlock + 2] = 1.0;
				if ( (x < 0.5) && (y < 0.5) && (z < 0.5) )
				{
					integralBlock1 += data[((i*testGrid[1]+j)*testGrid[0]+k)*nBlock + 0]*dV;
					integralBlock2 += data[((i*testGrid[1]+j)*testGrid[0]+k)*nBlock + 1]*dV;
					integralBlock3 += data[((i*testGrid[1]+j)*testGrid[0]+k)*nBlock + 2]*dV;
				}

				//sample the 1/2,1/2,1/2 part of the unit cell with a dense set of points
				//which form a regular mesh (irrelevant for the linear interpolation algorithm)
				testPoints[((i*testGrid[1]+j)*testGrid[0]+k)*3 + 0 ] = x/2.0;
				testPoints[((i*testGrid[1]+j)*testGrid[0]+k)*3 + 1 ] = y/2.0;
				testPoints[((i*testGrid[1]+j)*testGrid[0]+k)*3 + 2 ] = z/2.0;
			}

	//Check if our approximation results in a similar values as the analytic integral
	BOOST_REQUIRE( std::fabs(integralBlock1 - 1/std::pow(M_PI,3)) < accuracy );
	BOOST_REQUIRE( std::fabs(integralBlock2 -  1/std::pow(M_PI,3)) < accuracy );
	BOOST_REQUIRE( std::fabs(integralBlock3 - 1.0/8.0) < accuracy );

	elephon::LatticeStructure::RegularBareGrid regularGrid;
	regularGrid.initialize( testGrid );
	elephon::Algorithms::TrilinearInterpolation triLin(regularGrid);
	std::vector<int> indices;
	triLin.data_query(testPoints,indices);

	//copy the necessary indices
	std::vector<double> neededData(indices.size()*nBlock);
	for (int i = 0 ; i < indices.size(); ++i )
		for ( int ib = 0 ; ib < nBlock; ++ib)
			neededData[i*nBlock+ib] = data[indices[i]*nBlock+ib];

	std::vector<double> interpolatedData;
	triLin.interpolate(nBlock,neededData,interpolatedData);
	double dVInterpol = dV/8.0;
	double integral1 = 0.0;
	double integral2 = 0.0;
	double integral3 = 0.0;
	int nPoints = testPoints.size()/3;
	for ( int ik = 0 ; ik < nPoints; ++ik)
	{
		integral1 += interpolatedData[ik*nBlock + 0]*dVInterpol;
		integral2 += interpolatedData[ik*nBlock + 1]*dVInterpol;
		integral3 += interpolatedData[ik*nBlock + 2]*dVInterpol;
	}
	BOOST_REQUIRE( std::fabs(integral1 - integralBlock1) < accuracy);
	BOOST_REQUIRE( std::fabs(integral2 - integralBlock2) < accuracy);
	BOOST_REQUIRE( std::fabs(integral3 - integralBlock3) < accuracy);
}
