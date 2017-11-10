/*	This file test_ElectronPhononCoupling.cpp is part of elephon.
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
 *  Created on: Jul 4, 2017
 *      Author: A. Linscheid
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE PhononStructure
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "PhononStructure/ElectronPhononCoupling.h"
#include "ElectronicStructure/Wavefunctions.h"
#include "fixtures/MockStartup.h"
#include "fixtures/DataLoader.h"
#include "fixtures/scenarios.h"
#include <vector>

void
load_reference_data(boost::filesystem::path filename,
					std::vector<double> & qpoints,
					std::vector<int> & kdim,
					std::vector<double> & kpoints,
					int & numBands,
					int & numModes,
					std::vector<double> & gkkpData)
{
	std::ifstream file(filename.c_str());
	if ( ! file.good() )
		throw std::runtime_error(std::string("Problem opening file for reading: ")+filename.string());

	std::string line;
	std::getline(file, line); // comment

	{	// k sampling
		std::getline(file, line);
		std::stringstream ss(line);
		kdim.resize(3);
		ss >> kdim[0] >> kdim[1] >> kdim[2];
	}
	{	// q sampling
		std::getline(file, line);
		std::stringstream ss(line);
		qpoints.clear();
		double qi;
		while( ss >> qi )
		{
			qpoints.push_back(qi);
		}
	}
	{	// num bands
		std::getline(file, line);
		std::stringstream ss(line);
		ss >> numBands;
	}
	{	// num modes
		std::getline(file, line);
		std::stringstream ss(line);
		ss >> numModes;
	}

	int numK = kdim[0]*kdim[1]*kdim[2];
	kpoints.resize(numK*3);

	int numEle = qpoints.size()/3 * numK * numBands*numBands*numModes;
	gkkpData.resize(numEle);

	int c = 0, ck = 0;
	while(std::getline(file, line))
	{
		assert(c + numBands*numBands*numModes <= numEle);
		std::stringstream ss(line);
		if ( ck >= numK )
			ck = 0;
		for ( int i = 0 ; i < 3; ++i)
			ss >> kpoints[ck*3+i];
		for ( int ib = 0 ; ib < numBands*numBands*numModes; ++ib)
			ss >> gkkpData[c++];
	    ++ck;
	}
	assert(ck == numK);
	assert(c == numEle);

	const double mHaToEVEnergyConversion = 27.21138602/1000.0;
	for ( auto &g : gkkpData)
		g *= mHaToEVEnergyConversion*mHaToEVEnergyConversion;
}

BOOST_AUTO_TEST_CASE( Gkkp_generate_regular_k_grid_q_zero )
{
//	auto resHandl = elephon::test::fixtures::scenarios::load_Al_fcc_primitive_vasp_sc2x2x2();
//
//	// this reads in the reference data from the modified version of QE that writes these files
//	std::vector<double> qpoints;
//	std::vector<int> kdim;
//	std::vector<double> kpoints;
//	int numBands, numModes;
//	std::vector<double> gkkpDataRef;
//	auto ref_data_filename = boost::filesystem::path(resHandl->get_optns().get_root_dir())
//														/ "gkkp_reference.dat";
//	load_reference_data(ref_data_filename,
//						qpoints,
//						kdim,
//						kpoints,
//						numBands,
//						numModes,
//						gkkpDataRef);
//
//	std::vector<int> bandsList = { 1, 2 };
//	auto bandspList = bandsList;
//
//	int nq = qpoints.size()/3;
//	int nk = kpoints.size()/3;
//
//	elephon::PhononStructure::ElectronPhononCoupling gkkp;
//	double integral = 0.0, integralRef = 0.0;
//	std::vector<std::complex<float>> gkkpData;
//	for ( int iq = 0 ; iq < nq ; ++iq )
//	{
//		// generate_gkkp_and_phonon computes every combination of k, k'
//		// the reference data is for k0+q0, k0. Thus, we need to compute k'=k0, k=k0+q0
//		std::vector<double> kPrimePoints = kpoints;
//		for ( int ik = 0; ik < nk; ++ik)
//			for ( int i = 0; i < 3; ++i)
//			{
//				kPrimePoints[ik*3+i] += qpoints[iq*3+i];
//				kPrimePoints[ik*3+i] -= std::floor(kPrimePoints[ik*3+i]+0.5);
//			}
//
//		// PLEASE NOTE for manual inspection:
//		// The matrix elements computed by generate_gkkp_and_phonon are not symmetrized.
//		// Since the sum of the absolute value squared over the star of a given k point
//		// has the symmetry of the small group of q = k - k', they should be symmetrized
//		// by this small group. In the integral it does not matter, though.
//		gkkp.generate_gkkp_energy_units(
//				kPrimePoints, kpoints,
//				bandsList, bandspList,
//				resHandl->get_phonon_obj(),
//				resHandl->get_displacement_potential_obj(),
//				resHandl->get_wfct_obj(),
//				gkkpData);
//		BOOST_REQUIRE_EQUAL(gkkpData.size(), nk*numBands*numBands*numModes);
//
//		for ( int ik = 0 ; ik < nk ; ++ik )
//			for ( int ib = 0 ; ib < numBands ; ++ib)
//				for ( int ibp = 0 ; ibp < numBands ; ++ibp)
//					for ( int inu = 0 ; inu < numModes ; ++inu)
//					{
//						int cnsk = ((ik*numBands+ib)*numBands+ibp)*numModes + inu;
//						int cnsq = (((iq*nk+ik)*numBands+ib)*numBands+ibp)*numModes + inu;
//						integral += std::real(std::conj(gkkpData[cnsk])*(gkkpData[cnsk]));
//						integralRef += gkkpDataRef[cnsq] ;
//					}
//	}
//
//	BOOST_CHECK_CLOSE(integral, integralRef, 0.001);
}

