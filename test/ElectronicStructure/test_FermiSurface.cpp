/*	This file test_FermiSurface.cpp is part of elephon.
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
 *  Created on: Apr 27, 2017
 *      Author: A. Linscheid
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Input_test
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "ElectronicStructure/FermiSurface.h"
#include "ElectronicStructure/GradientFFTReciprocalGrid.h"
#include "Algorithms/TrilinearInterpolation.h"
#include "LatticeStructure/RegularBareGrid.h"
#include "fixtures/MockStartup.h"
#include <assert.h>
#include <vector>
#include <fstream>
#include <exception>
#include <algorithm>
#include <iostream>

/**
 * This test checks if the surface area of a sphere is correctly approximated
 */
BOOST_AUTO_TEST_CASE( Check_surface_sphere )
{
	//construct data which will yield a sphere of radius r
	double const r = 0.25;
	elephon::LatticeStructure::RegularBareGrid grid;
	grid.initialize( {100,100,100}, true, 1e-6, {0.0, 0.0, 0.0}, elephon::LatticeStructure::LatticeModule() );
	auto d = grid.get_grid_dim();
	std::vector<double> data(d[0]*d[1]*d[2]);
	for (int k = 0 ; k < d[2]; ++k)
	  for (int j = 0 ; j < d[1]; ++j)
			for (int i = 0 ; i < d[0]; ++i)
		  {
			  double x = (i<d[0]/2?double(i):double(i)-d[0])/double(d[0]);
			  double y = (j<d[1]/2?double(j):double(j)-d[1])/double(d[1]);
			  double z = (k<d[2]/2?double(k):double(k)-d[2])/double(d[2]);
			  data[(k*d[1]+j)*d[0]+i] = std::sqrt(x*x+y*y+z*z)/r;
		  }


	int targetNumPoints = 10000;
	elephon::ElectronicStructure::FermiSurface fs;
	fs.triangulate_surface(
			grid,
			1,
			data,
			targetNumPoints,
			1.0);

	double surfaceArea = 0;

	std::vector<double> kfweights = fs.get_Fermi_weights_for_band(0);
	for (int ikf = 0 ; ikf < kfweights.size(); ++ikf)
		surfaceArea += kfweights[ikf]/std::pow(2*M_PI,2);

	std::cout << "Surface area of a sphere (triangulated): "
			<< surfaceArea << " expected: " << 4*M_PI*r*r << std::endl;
	BOOST_REQUIRE( std::fabs(surfaceArea-4*M_PI*r*r) < 1e-2);
}

/**
 * This test compares the numerical results for the DOS against the analytically
 * solvable free electron gas.
 */
BOOST_AUTO_TEST_CASE( free_electron_gas )
{
	std::cout << "Compare numercial results against analytical results for the free e- gas " << std::endl;
	//These are the parameters that define the free electron gas model
	double const meStar = 0.200;
	double const Ne = 0.001;

	double const Ef = 1.0/(2.0*meStar)*std::pow(3*M_PI*M_PI*Ne,2.0/3.0);
	elephon::LatticeStructure::RegularBareGrid grid;
	grid.initialize( {100,100,100}, true, 1e-6, {0.0, 0.0, 0.0}, elephon::LatticeStructure::LatticeModule() );
	std::vector<int> d = grid.get_grid_dim();
	std::vector<double> data(d[0]*d[1]*d[2]);
	for (int k = 0 ; k < d[2]; ++k)
		for (int j = 0 ; j < d[1]; ++j)
			for (int i = 0 ; i < d[0]; ++i)
			  {
				  double kx = (i<d[0]/2?double(i):double(i)-d[0])/double(d[0]);
				  double ky = (j<d[1]/2?double(j):double(j)-d[1])/double(d[1]);
				  double kz = (k<d[2]/2?double(k):double(k)-d[2])/double(d[2]);
				  data[(k*d[1]+j)*d[0]+i] = (kx*kx+ky*ky+kz*kz)/(2.0*meStar)-Ef;
			  }

	elephon::ElectronicStructure::GradientFFTReciprocalGrid gradE;
	gradE.compute_gradient(
			d,
			grid.get_lattice(),
			1,
			data);

	elephon::Algorithms::TrilinearInterpolation triLin(grid);

	//Here we compute if the DOS is computed correctly for the analytically solvable model
	const int numESteps = 20;
	const double estart = *std::min_element(data.begin(),data.end());//eV
	const double eend = 0.0;//eV
	const double de =  (eend-estart)/double(numESteps);
	for ( int ie = 0 ; ie < numESteps; ++ie)
	{
		double e = estart + de*ie;

		int targetNumPoints = 1000;
		elephon::ElectronicStructure::FermiSurface fs;
		fs.triangulate_surface(
				grid,
				1,
				data,
				targetNumPoints,
				e);

		//Interpolate the velocities to the Fermi level
		std::vector<int> reqestedIndices;
		std::vector<double> kfVect = fs.get_Fermi_vectors_for_band(0);

		if ( kfVect.size() == 0 )
			continue;

		triLin.data_query( kfVect, reqestedIndices );

		std::vector<double> gradDataAtRequestedIndices;
		gradE.copy_data( reqestedIndices, std::vector<int>({0}), gradDataAtRequestedIndices );

		std::vector<double> FermiVelocities;
		triLin.interpolate(3,gradDataAtRequestedIndices,FermiVelocities);

		BOOST_REQUIRE( FermiVelocities.size() == kfVect.size() );

		std::vector<double> kfweights = fs.get_Fermi_weights_for_band(0);
		double dos = 0;
		for (int ikf = 0 ; ikf < kfVect.size()/3; ++ikf)
		{
			double modGradE = std::sqrt(std::pow(FermiVelocities[ikf*3+0],2)
									+std::pow(FermiVelocities[ikf*3+1],2)
									+std::pow(FermiVelocities[ikf*3+2],2));
			dos += kfweights[ikf]/modGradE/std::pow(2*M_PI,2);
		}

		//per spin
		double analytic_dos = 1.0/std::pow(2*M_PI,2)*std::pow(2*meStar,1.5)*std::sqrt(e+Ef);
		BOOST_REQUIRE( (std::fabs(dos-analytic_dos)/(dos+analytic_dos) < 1e-2)
				or (analytic_dos < 1e-1) );
	}
}

/**
 * This test checks if the surface area of a sphere plus a spheroid is correctly approximated
 */
BOOST_AUTO_TEST_CASE( Check_surface_sphere_plus_spheroid )
{
	std::cout << "Compute the surface area of a sphere plus spheroid " << std::endl;
	//construct data which will yield a sphere of radius r
	double const r1 = 0.25;
	double const a = 0.25;
	double const c = 0.35;
	const int nBnd = 2;
	elephon::LatticeStructure::RegularBareGrid grid;
	grid.initialize( {100,100,100}, true, 1e-6, {0.0, 0.0, 0.0}, elephon::LatticeStructure::LatticeModule() );
	std::vector<int> d = grid.get_grid_dim();
	std::vector<double> data(d[0]*d[1]*d[2]*nBnd);
	for (int k = 0 ; k < d[2]; ++k)
		for (int j = 0 ; j < d[1]; ++j)
			for (int i = 0 ; i < d[0]; ++i)
			{
			  double x = (i<d[0]/2?double(i):double(i)-d[0])/double(d[0]);
			  double y = (j<d[1]/2?double(j):double(j)-d[1])/double(d[1]);
			  double z = (k<d[2]/2?double(k):double(k)-d[2])/double(d[2]);
			  data[((k*d[1]+j)*d[0]+i)*2] = (x*x+y*y+z*z)/(r1*r1);
			  data[((k*d[1]+j)*d[0]+i)*2+1] = (x*x+y*y)/(a*a)+z*z/(c*c);
			}

	elephon::LatticeStructure::LatticeModule lattice;

	int targetNumPoints = 10000;
	elephon::ElectronicStructure::FermiSurface fs;
	fs.triangulate_surface(
			grid,
			nBnd,
			data,
			targetNumPoints,
			1.0);

	double surfaceArea = 0;
	for ( int ib = 0 ; ib < nBnd; ++ib)
	{
		std::vector<double> kfweights = fs.get_Fermi_weights_for_band(ib);
		for (int ikf = 0 ; ikf < kfweights.size(); ++ikf)
			surfaceArea += kfweights[ikf]/std::pow(2*M_PI,2);
	}

	const double sphereArea = 4*M_PI*r1*r1;
	const double e = std::sqrt(1-a*a/(c*c));
	const double spheroidArea = 2*M_PI*a*a*(1+c/(e*a)*std::asin(e));

	BOOST_REQUIRE( std::fabs(surfaceArea-spheroidArea-sphereArea) < 2e-3);
}

/**
 * This test computes the DOS from a straight forward delta function approximation
 * and compares the (better) DOS from the energy surface trianglulation against it.
 * Test uses a two band, 2D cosine band structure (next nearest neighbour hopping)
 */
BOOST_AUTO_TEST_CASE( Check_DOS_2BdndCos )
{
	std::cout << "Compute DOS of a cosine band model " << std::endl;
	//These are the parameters that define the cosine model
	double const W1 = 200; // electron
	double const W2 = 100; // hole
	double const Ee = -50;
	double const Eg = 10;
	int nBnd = 2;

	elephon::LatticeStructure::RegularBareGrid grid;
	grid.initialize(  {99,101,100}, true, 1e-6, {0.0, 0.0, 0.0}, elephon::LatticeStructure::LatticeModule() );
	std::vector<int> d = grid.get_grid_dim();
	std::vector<double> data(d[0]*d[1]*d[2]*nBnd);
	for (int k = 0 ; k < d[2]; ++k)
		for (int j = 0 ; j < d[1]; ++j)
			for (int i = 0 ; i < d[0]; ++i)
			{
				double kx = (i<d[0]/2?double(i):double(i)-d[0])/double(d[0]);
				double ky = (j<d[1]/2?double(j):double(j)-d[1])/double(d[1]);
				data[((k*d[1]+j)*d[0]+i)*nBnd+0] =
					  W1*(std::cos( 2*M_PI*kx )+std::cos( 2*M_PI*ky )-2.0)+Ee;
				data[((k*d[1]+j)*d[0]+i)*nBnd+1] =
					  W2*(std::cos( 2*M_PI*kx )+std::cos( 2*M_PI*ky )-2.0)+Eg;
			}

	//here we set parameter for the DOS caluation that we cross-check
	const int numESteps = 50;
	const double estart = *std::min_element(data.begin(),data.end());//eV
	const double eend = *std::max_element(data.begin(),data.end());//eV
	const double de =  (eend-estart)/double(numESteps);

	//Here we use a very simple delta function
	//approximation to estimate the DOS
	std::vector<double> dos_ref(numESteps, 0.0 );
	double dKV = 1.0/double(d[0]*d[1]*d[2]);//Unit volume
	for (int i = 0 ; i < numESteps; ++i)
	{
		double eWinMin = estart+de*(double(i)-1.0);
		double eWinMax = estart+de*(double(i)+0.0);
		for ( auto e : data )
			if ( e >= eWinMin and e < eWinMax )
				dos_ref[i] += dKV/de;
	}

	elephon::ElectronicStructure::GradientFFTReciprocalGrid gradE;
	gradE.compute_gradient(
			d,
			grid.get_lattice(),
			nBnd,
			data);

	elephon::Algorithms::TrilinearInterpolation triLin(grid);

	double Nelec = 0;

	//Here we check if the DOS is computed correctly against the simple reference
	for ( int ie = 0 ; ie < numESteps; ++ie)
	{
		double e = estart + de*ie;

		int targetNumPoints = 10000;
		elephon::ElectronicStructure::FermiSurface fs;
		fs.triangulate_surface(
				grid,
				nBnd,
				data,
				targetNumPoints,
				e);

		double dos = 0;
		for ( int ib = 0 ; ib < nBnd; ++ib)
		{
			//Interpolate the velocities to the Fermi level
			std::vector<int> reqestedIndices;
			std::vector<double> kfVect = fs.get_Fermi_vectors_for_band(ib);

			if ( kfVect.size() == 0 )
				continue;

			triLin.data_query( kfVect, reqestedIndices );

			std::vector<double> gradDataAtRequestedIndices;
			gradE.copy_data( reqestedIndices, std::vector<int>({ib}), gradDataAtRequestedIndices );

			std::vector<double> FermiVelocities;
			triLin.interpolate(3,gradDataAtRequestedIndices,FermiVelocities);

			BOOST_REQUIRE( FermiVelocities.size() == kfVect.size() );

			std::vector<double> kfweights = fs.get_Fermi_weights_for_band(ib);
			for (int ikf = 0 ; ikf < kfVect.size()/3; ++ikf)
			{
				double modGradE = std::sqrt(std::pow(FermiVelocities[ikf*3+0],2)
										+std::pow(FermiVelocities[ikf*3+1],2)
										+std::pow(FermiVelocities[ikf*3+2],2));
				dos += kfweights[ikf]/modGradE/std::pow(2*M_PI,3);// by a where a is 1
			}
		}
		if ( e <= 0 )
			Nelec += de*dos;
	}

	double Nelec_ref = 0;
	for ( int i = 0 ; i < numESteps; ++i )
		if ( estart + de*i <= 0 )
			Nelec_ref += de*dos_ref[i];

	//Here we check the integral
	std::cout << "Almost filled 2-band cosine model electrons are (test,ref): " << Nelec
			<<'\t'<< Nelec_ref << std::endl;
	BOOST_REQUIRE( std::fabs(Nelec-Nelec_ref) <= 1e-1 );
}

/**
 * This test computes the DOS from a straight forward delta function approximation
 * and compares the (better) DOS from the energy surface triangulation against it.
 * This check loads a data file for LiFeAs and computes the LiFeAs DOS.
 */
BOOST_AUTO_TEST_CASE( Check_DOS_LiFeAs )
{
	std::cout << "Compute the DOS of a LiFeAs with 10 bands " << std::endl;
	//load data from the LiFeAs electronic structure
	test::fixtures::MockStartup ms;
	auto testd = ms.get_data_for_testing_dir();
	std::ifstream input( (testd / "LiFeAs_energies.dat").c_str(), std::ios::in | std::ios::binary );

	if (not input)
		throw std::runtime_error ( std::string("Hard coded file ")+(testd / "LiFeAs_energies.dat").string()+" not readable ");

	std::streampos fileSize;
	int sizeOfBuffer;

	// Get the size of the file
	input.seekg(0, std::ios::end);
	fileSize = input.tellg();
	input.seekg(0, std::ios::beg);

	sizeOfBuffer = fileSize / sizeof(double);
	std::vector<double> buffer(sizeOfBuffer);

	input.read(reinterpret_cast<char*>(buffer.data()), fileSize);

	int nBnd = std::round(buffer[0]);
	int nkx = std::round(buffer[1]);
	int nky = std::round(buffer[2]);
	int nkz = std::round(buffer[3]);

	elephon::LatticeStructure::LatticeModule lattice(
										{   7.050131 , 0.000000 , 0.000000 ,
											0.000000 , 7.050131 , 0.000000 ,
											0.000000 , 0.000000 , 11.45919 } );

	elephon::LatticeStructure::RegularBareGrid grid;
	grid.initialize({nkx,nky,nkz}, true, 1e-6, {0.0, 0.0, 0.0}, lattice );
	std::vector<int> d = grid.get_grid_dim();
	assert( nBnd*nkx*nky*nkz == (buffer.size()-4) );

	std::vector<double> energies(&buffer[4],&buffer[4]+nBnd*nkx*nky*nkz);


	elephon::ElectronicStructure::GradientFFTReciprocalGrid gradE;
	gradE.compute_gradient(d,
			lattice,
			nBnd,
			energies);

	//Compute DOS
	const int numESteps = 100;

	const double estart = *std::min_element(energies.begin(),energies.end());//eV
	const double eend = 0;//eV
	std::vector<double> dos(numESteps,0.0);
	double de =  (eend-estart)/double(numESteps);

	//Here we use a very simple delta function
	//approximation to estimate the DOS
	std::vector<double> dos_ref(numESteps, 0.0 );
	for (int i = 0 ; i < numESteps; ++i)
	{
		double eWinMin = estart+de*(double(i)-0.5);
		double eWinMax = estart+de*(double(i)+0.5);
		for ( auto e : energies )
			if ( e >= eWinMin and e < eWinMax )
				dos_ref[i] += 1.0/de/double(d[0]*d[1]*d[2]);
	}

	elephon::Algorithms::TrilinearInterpolation triLin(grid);

	for ( int ie = 0 ; ie < numESteps; ++ie)
	{
		double e = estart + de*ie;

		int targetNumPoints = 10000;
		elephon::ElectronicStructure::FermiSurface fs;
		fs.triangulate_surface(
				grid,
				nBnd,
				energies,
				targetNumPoints,
				e);

		for ( int ib = 0 ; ib < nBnd; ++ib)
		{
			//Interpolate the velocities to the Fermi level
			std::vector<int> reqestedIndices;
			std::vector<double> kfVect = fs.get_Fermi_vectors_for_band(ib);

			if ( kfVect.size() == 0 )
				continue;

			triLin.data_query( kfVect, reqestedIndices );;

			std::vector<double> gradDataAtRequestedIndices;
			gradE.copy_data( reqestedIndices, std::vector<int>({ib}), gradDataAtRequestedIndices );

			std::vector<double> FermiVelocities;
			triLin.interpolate(3,gradDataAtRequestedIndices,FermiVelocities);

			BOOST_REQUIRE( FermiVelocities.size() == kfVect.size() );

			std::vector<double> kfweights = fs.get_Fermi_weights_for_band(ib);
			for (int ikf = 0 ; ikf < kfVect.size()/3; ++ikf)
			{
				double modGradE = std::sqrt(std::pow(FermiVelocities[ikf*3+0],2)
										+std::pow(FermiVelocities[ikf*3+1],2)
										+std::pow(FermiVelocities[ikf*3+2],2));
				dos[ie] += kfweights[ikf]/modGradE*lattice.get_volume()/std::pow(2*M_PI,3);
			}
		}
	}

	//Check the number of electrons in the system
	double Nelec = 0;
	double Nelec_ref = 0;
	for ( int i = 0 ; i < numESteps; ++i )
		if ( estart + de*i <= 0 )
		{
			Nelec += de*dos[i];
			Nelec_ref += de*dos_ref[i];
		}
	std::cout << "Computed number of electrons (simple, surface-integral method): "
			<< Nelec_ref << '\t' << Nelec << std::endl;

	BOOST_REQUIRE( std::fabs(Nelec-Nelec_ref)/(Nelec+Nelec_ref) < 1e-1);
}
