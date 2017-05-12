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
	std::vector<size_t> grid = {100,100,100};
	std::vector<double> data(grid[0]*grid[1]*grid[2]);
	for (size_t i = 0 ; i < grid[0]; ++i)
	  for (size_t j = 0 ; j < grid[1]; ++j)
		  for (size_t k = 0 ; k < grid[2]; ++k)
		  {
			  double x = (i<grid[0]/2?double(i):double(i)-grid[0])/double(grid[0]);
			  double y = (j<grid[1]/2?double(j):double(j)-grid[1])/double(grid[1]);
			  double z = (k<grid[2]/2?double(k):double(k)-grid[2])/double(grid[2]);
			  data[(i*grid[1]+j)*grid[2]+k] = std::sqrt(x*x+y*y+z*z)/r;
		  }

	std::vector<double> reciprocalLatticeMatrix
											= {   1.000000 , 0.000000 , 0.000000 ,
												  0.000000 , 1.000000 , 0.000000 ,
												  0.000000 , 0.000000 , 1.000000 	};

	size_t targetNumPoints = 10000;
	elephon::ElectronicStructure::FermiSurface fs;
	fs.triangulate_surface(
			grid,
			reciprocalLatticeMatrix,
			1,
			data,
			targetNumPoints,
			1.0);

	double surfaceArea = 0;

	std::vector<double> kfweights = fs.get_Fermi_weights_for_band(0);
	for (size_t ikf = 0 ; ikf < kfweights.size(); ++ikf)
		surfaceArea += kfweights[ikf];

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
	//These are the parameters that define the free electron gas model
	double const meStar = 0.200;
	double const Ne = 0.001;

	double const Ef = 1.0/(2.0*meStar)*std::pow(3*M_PI*M_PI*Ne,2.0/3.0);
	std::vector<size_t> grid = {100,100,100};
	std::vector<double> data(grid[0]*grid[1]*grid[2]);
	for (size_t i = 0 ; i < grid[0]; ++i)
	  for (size_t j = 0 ; j < grid[1]; ++j)
		  for (size_t k = 0 ; k < grid[2]; ++k)
		  {
			  double kx = (i<grid[0]/2?double(i):double(i)-grid[0])/double(grid[0]);
			  double ky = (j<grid[1]/2?double(j):double(j)-grid[1])/double(grid[1]);
			  double kz = (k<grid[2]/2?double(k):double(k)-grid[2])/double(grid[2]);
			  data[(i*grid[1]+j)*grid[2]+k] = (kx*kx+ky*ky+kz*kz)/(2.0*meStar)-Ef;
		  }

	std::vector<double> latticeMatrix = {   1.000000 , 0.000000 , 0.000000 ,
											0.000000 , 1.000000 , 0.000000 ,
											0.000000 , 0.000000 , 1.000000 };

	//for the analytic model the k vectors is measured in units of 1, not 2pi
	std::vector<double> reciprocalLatticeMatrix
									  = {   1.000000 , 0.000000 , 0.000000 ,
											0.000000 , 1.000000 , 0.000000 ,
											0.000000 , 0.000000 , 1.000000 };

	elephon::ElectronicStructure::GradientFFTReciprocalGrid gradE;
	gradE.compute_gradient(grid,
			latticeMatrix,
			1,
			data);

	//Here we compute if the DOS is computed correctly for the analytically solvable model
	const size_t numESteps = 20;
	const double estart = *std::min_element(data.begin(),data.end());//eV
	const double eend = 0.0;//eV
	const double de =  (eend-estart)/double(numESteps);
	for ( size_t ie = 0 ; ie < numESteps; ++ie)
	{
		double e = estart + de*ie;

		size_t targetNumPoints = 10000;
		elephon::ElectronicStructure::FermiSurface fs;
		fs.triangulate_surface(
				grid,
				reciprocalLatticeMatrix,
				1,
				data,
				targetNumPoints,
				e);

		//Interpolate the velocities to the Fermi level
		elephon::Algorithms::TrilinearInterpolation triLin(grid);
		std::vector<size_t> reqestedIndices;
		std::vector<double> kfVect = fs.get_Fermi_vectors_for_band(0);

		if ( kfVect.size() == 0 )
			continue;

		triLin.data_query( kfVect, reqestedIndices );

		std::vector<double> gradDataAtRequestedIndices;
		gradE.copy_data( reqestedIndices, std::vector<size_t>({0}), gradDataAtRequestedIndices );

		std::vector<double> FermiVelocities;
		triLin.interpolate(3,gradDataAtRequestedIndices,FermiVelocities);

		BOOST_REQUIRE( FermiVelocities.size() == kfVect.size() );

		std::vector<double> kfweights = fs.get_Fermi_weights_for_band(0);
		double dos = 0;
		for (size_t ikf = 0 ; ikf < kfVect.size()/3; ++ikf)
		{
			double modGradE = std::sqrt(std::pow(FermiVelocities[ikf*3+0],2)
									+std::pow(FermiVelocities[ikf*3+1],2)
									+std::pow(FermiVelocities[ikf*3+2],2));
			dos += kfweights[ikf]/modGradE;
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
	//construct data which will yield a sphere of radius r
	double const r1 = 0.25;
	double const a = 0.25;
	double const c = 0.35;
	const size_t nBnd = 2;
	std::vector<size_t> grid = {100,100,100};
	std::vector<double> data(grid[0]*grid[1]*grid[2]*nBnd);
	for (size_t i = 0 ; i < grid[0]; ++i)
	  for (size_t j = 0 ; j < grid[1]; ++j)
		  for (size_t k = 0 ; k < grid[2]; ++k)
		  {
			  double x = (i<grid[0]/2?double(i):double(i)-grid[0])/double(grid[0]);
			  double y = (j<grid[1]/2?double(j):double(j)-grid[1])/double(grid[1]);
			  double z = (k<grid[2]/2?double(k):double(k)-grid[2])/double(grid[2]);
			  data[((i*grid[1]+j)*grid[2]+k)*2] = (x*x+y*y+z*z)/(r1*r1);
			  data[((i*grid[1]+j)*grid[2]+k)*2+1] = (x*x+y*y)/(a*a)+z*z/(c*c);
		  }

	std::vector<double> reciprocalLatticeMatrix
											= {   1.000000 , 0.000000 , 0.000000 ,
												  0.000000 , 1.000000 , 0.000000 ,
												  0.000000 , 0.000000 , 1.000000 	};

	size_t targetNumPoints = 10000;
	elephon::ElectronicStructure::FermiSurface fs;
	fs.triangulate_surface(
			grid,
			reciprocalLatticeMatrix,
			nBnd,
			data,
			targetNumPoints,
			1.0);

	double surfaceArea = 0;
	for ( size_t ib = 0 ; ib < nBnd; ++ib)
	{
		std::vector<double> kfweights = fs.get_Fermi_weights_for_band(ib);
		for (size_t ikf = 0 ; ikf < kfweights.size(); ++ikf)
			surfaceArea += kfweights[ikf];
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
	//These are the parameters that define the cosine model
	double const W1 = 200; // electron
	double const W2 = 100; // hole
	double const Ee = -50;
	double const Eg = 10;
	size_t nBnd = 2;

	std::vector<size_t> grid = {200,201,100};
	std::vector<double> data(grid[0]*grid[1]*grid[2]*nBnd);
	for (size_t i = 0 ; i < grid[0]; ++i)
	  for (size_t j = 0 ; j < grid[1]; ++j)
		  for (size_t k = 0 ; k < grid[2]; ++k)
		  {
			  double kx = (i<grid[0]/2?double(i):double(i)-grid[0])/double(grid[0]);
			  double ky = (j<grid[1]/2?double(j):double(j)-grid[1])/double(grid[1]);
			  data[((i*grid[1]+j)*grid[2]+k)*nBnd+0] =
					  W1*(std::cos( 2*M_PI*kx )+std::cos( 2*M_PI*ky )-2.0)+Ee;
			  data[((i*grid[1]+j)*grid[2]+k)*nBnd+1] =
					  W2*(std::cos( 2*M_PI*kx )+std::cos( 2*M_PI*ky )-2.0)+Eg;
		  }

	//here we set parameter for the DOS caluation that we cross-check
	const size_t numESteps = 20;
	const double estart = *std::min_element(data.begin(),data.end());//eV
	const double eend = *std::max_element(data.begin(),data.end());//eV
	const double de =  (eend-estart)/double(numESteps);

	//Here we use a very simple delta function
	//approximation to estimate the DOS
	std::vector<double> dos_ref(numESteps, 0.0 );
	double dKV = 1.0/double(grid[0]*grid[1]*grid[2]);//Unit volume
	for (size_t i = 0 ; i < numESteps; ++i)
	{
		double eWinMin = estart+de*(double(i)-1.0);
		double eWinMax = estart+de*(double(i)+0.0);
		for ( auto e : data )
			if ( e >= eWinMin and e < eWinMax )
				dos_ref[i] += dKV/de;
	}

	std::vector<double> latticeMatrix = {   1.000000 , 0.000000 , 0.000000 ,
											0.000000 , 1.000000 , 0.000000 ,
											0.000000 , 0.000000 , 1.000000 };
	std::vector<double> reciprocalLatticeMatrix
									  = {   2.0*M_PI , 0.000000 , 0.000000 ,
											0.000000 , 2.0*M_PI , 0.000000 ,
											0.000000 , 0.000000 , 2.0*M_PI };

	elephon::ElectronicStructure::GradientFFTReciprocalGrid gradE;
	gradE.compute_gradient(grid,
			latticeMatrix,
			nBnd,
			data);

	double Nelec = 0;

	//Here we check if the DOS is computed correctly against the simple reference
	for ( size_t ie = 0 ; ie < numESteps; ++ie)
	{
		double e = estart + de*ie;

		size_t targetNumPoints = 10000;
		elephon::ElectronicStructure::FermiSurface fs;
		fs.triangulate_surface(
				grid,
				reciprocalLatticeMatrix,
				nBnd,
				data,
				targetNumPoints,
				e);

		double dos = 0;
		for ( size_t ib = 0 ; ib < nBnd; ++ib)
		{
			//Interpolate the velocities to the Fermi level
			elephon::Algorithms::TrilinearInterpolation triLin(grid);
			std::vector<size_t> reqestedIndices;
			std::vector<double> kfVect = fs.get_Fermi_vectors_for_band(ib);

			if ( kfVect.size() == 0 )
				continue;

			triLin.data_query( kfVect, reqestedIndices );

			std::vector<double> gradDataAtRequestedIndices;
			gradE.copy_data( reqestedIndices, std::vector<size_t>({ib}), gradDataAtRequestedIndices );

			std::vector<double> FermiVelocities;
			triLin.interpolate(3,gradDataAtRequestedIndices,FermiVelocities);

			BOOST_REQUIRE( FermiVelocities.size() == kfVect.size() );

			std::vector<double> kfweights = fs.get_Fermi_weights_for_band(ib);
			for (size_t ikf = 0 ; ikf < kfVect.size()/3; ++ikf)
			{
				double modGradE = std::sqrt(std::pow(FermiVelocities[ikf*3+0],2)
										+std::pow(FermiVelocities[ikf*3+1],2)
										+std::pow(FermiVelocities[ikf*3+2],2));
				dos += kfweights[ikf]/modGradE/std::pow(2*M_PI,3);
			}
		}
		if ( e <= 0 )
			Nelec += de*dos;
	}

	double Nelec_ref = 0;
	for ( size_t i = 0 ; i < numESteps; ++i )
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
	//load data from the LiFeAs electronic structure
	boost::filesystem::path p(__FILE__);
	boost::filesystem::path dir = p.parent_path();
	std::string LiFeAs_elstr_f = std::string(dir.c_str())+"/../ElectronicStructure/LiFeAs_energies.dat";
	std::ifstream input( LiFeAs_elstr_f.c_str(), std::ios::in | std::ios::binary );

	if (not input)
		throw std::logic_error ( std::string("Hard coded file ")+LiFeAs_elstr_f+" not readable ");

	std::streampos fileSize;
	size_t sizeOfBuffer;

	// Get the size of the file
	input.seekg(0, std::ios::end);
	fileSize = input.tellg();
	input.seekg(0, std::ios::beg);

	sizeOfBuffer = fileSize / sizeof(double);
	std::vector<double> buffer(sizeOfBuffer);

	input.read(reinterpret_cast<char*>(buffer.data()), fileSize);

	size_t nBnd = std::round(buffer[0]);
	size_t nkx = std::round(buffer[1]);
	size_t nky = std::round(buffer[2]);
	size_t nkz = std::round(buffer[3]);
	std::vector<size_t> grid({nkx,nky,nkz});
	assert( nBnd*nkx*nky*nkz == (buffer.size()-4) );

	std::vector<double> energies(nBnd*nkx*nky*nkz);
	for (size_t i = 0 ; i < grid[0]; ++i)
		for (size_t j = 0 ; j < grid[1]; ++j)
			for (size_t k = 0 ; k < grid[2]; ++k)
				for (size_t ib = 0 ; ib < nBnd; ++ib)
				{
					size_t consqFortran = ((k*grid[1]+j)*grid[0] +i)*nBnd+ib;
					size_t consqCpp =  ((i* grid[1]+j)*grid[2]+k)*nBnd+ib;
					energies[consqCpp] = buffer[4+consqFortran];
				}

	std::vector<double> latticeMatrix = {   7.050131 , 0.000000 , 0.000000 ,
											0.000000 , 7.050131 , 0.000000 ,
											0.000000 , 0.000000 , 11.45919 };

	//units are 2pi/a
	std::vector<double> reciprocalLatticeMatrix
									  = {   0.891215 , 0.000000 , 0.000000 ,
											0.000000 , 0.891215 , 0.000000 ,
											0.000000 , 0.000000 , 0.548310 };

	double Volume = latticeMatrix[0]*latticeMatrix[4]*latticeMatrix[8];

	elephon::ElectronicStructure::GradientFFTReciprocalGrid gradE;
	gradE.compute_gradient(grid,
			latticeMatrix,
			nBnd,
			energies);

	//Compute DOS
	const size_t numESteps = 50;

	const double estart = *std::min_element(energies.begin(),energies.end());//eV
	const double eend = 0;//eV
	std::vector<double> dos(numESteps,0.0);
	double de =  (eend-estart)/double(numESteps);

	//Here we use a very simple delta function
	//approximation to estimate the DOS
	std::vector<double> dos_ref(numESteps, 0.0 );
	for (size_t i = 0 ; i < numESteps; ++i)
	{
		double eWinMin = estart+de*(double(i)-0.5);
		double eWinMax = estart+de*(double(i)+0.5);
		for ( auto e : energies )
			if ( e >= eWinMin and e < eWinMax )
				dos_ref[i] += 1.0/de/double(grid[0]*grid[1]*grid[2]);
	}

	for ( size_t ie = 0 ; ie < numESteps; ++ie)
	{
		double e = estart + de*ie;

		size_t targetNumPoints = 10000;
		elephon::ElectronicStructure::FermiSurface fs;
		fs.triangulate_surface(
				grid,
				reciprocalLatticeMatrix,
				nBnd,
				energies,
				targetNumPoints,
				e);


		for ( size_t ib = 0 ; ib < nBnd; ++ib)
		{
			//Interpolate the velocities to the Fermi level
			elephon::Algorithms::TrilinearInterpolation triLin(grid);
			std::vector<size_t> reqestedIndices;
			std::vector<double> kfVect = fs.get_Fermi_vectors_for_band(ib);

			if ( kfVect.size() == 0 )
				continue;

			triLin.data_query( kfVect, reqestedIndices );;

			std::vector<double> gradDataAtRequestedIndices;
			gradE.copy_data( reqestedIndices, std::vector<size_t>({ib}), gradDataAtRequestedIndices );

			std::vector<double> FermiVelocities;
			triLin.interpolate(3,gradDataAtRequestedIndices,FermiVelocities);

			BOOST_REQUIRE( FermiVelocities.size() == kfVect.size() );

			std::vector<double> kfweights = fs.get_Fermi_weights_for_band(ib);
			for (size_t ikf = 0 ; ikf < kfVect.size()/3; ++ikf)
			{
				double modGradE = std::sqrt(std::pow(FermiVelocities[ikf*3+0],2)
										+std::pow(FermiVelocities[ikf*3+1],2)
										+std::pow(FermiVelocities[ikf*3+2],2));
				dos[ie] += kfweights[ikf]/modGradE*Volume/std::pow(2*M_PI,3);
			}
		}
	}

	//Check the number of electrons in the system
	double Nelec = 0;
	double Nelec_ref = 0;
	for ( size_t i = 0 ; i < numESteps; ++i )
		if ( estart + de*i <= 0 )
		{
			Nelec += de*dos[i];
			Nelec_ref += de*dos_ref[i];
		}
	std::cout << "Computed number of electrons (simple, surface-integral method): "
			<< Nelec_ref << '\t' << Nelec << std::endl;

	BOOST_REQUIRE( std::fabs(Nelec-Nelec_ref)/(Nelec+Nelec_ref) < 1e-1);
}
