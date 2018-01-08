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
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "PhononStructure/ElectronPhononCoupling.h"
#include "ElectronicStructure/Wavefunctions.h"
#include "fixtures/MockStartup.h"
#include "fixtures/DataLoader.h"
#include "fixtures/scenarios.h"
#include <vector>
#include <math.h>
#include <algorithm>

BOOST_AUTO_TEST_SUITE( ElectronPhononCoupling )

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

	for ( auto &ki : kpoints)
		ki -= std::floor(ki +0.5);

	const double RyToEVEnergyConversion = 27.21138602/2.0;
	for ( auto &g : gkkpData)
		g *= RyToEVEnergyConversion*RyToEVEnergyConversion;
}

//Again, we have to outcomment this cross check below, because the time a test takes is unaccatable in debug mode.
BOOST_AUTO_TEST_CASE( Gkkp_generate_regular_k_grid_q_zero )
{
//	auto resHandl = elephon::test::fixtures::scenarios::load_Al_fcc_primitive_vasp_sc4x4x4();
//
//	// this reads in the reference data from the modified version of QE that writes these files
//	std::vector<double> qpoints;
//	std::vector<int> kdim;
//	std::vector<double> kpoints;
//	int numBands, numModes;
//	std::vector<double> gkkpDataRef;
//	auto ref_data_filename = boost::filesystem::path(resHandl->get_optns().get_root_dir())
//														/ "gkkp_reference.dat";
//	auto kgrid = resHandl->get_wfct_obj()->get_k_grid().view_bare_grid();
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
//	std::vector<int> kpointIndices;
//	std::vector<int> qpointIndex;
//	kgrid.get_list_reducible_lattice_point_indices(qpoints, qpointIndex);
//	kgrid.get_list_reducible_lattice_point_indices(kpoints, kpointIndices);
//	std::vector<std::pair<int,int>> kandKPPointIndices(kpointIndices.size());
//	std::vector<int> xyz(3), xyzQ(3);
//	kgrid.get_reducible_to_xyz(qpointIndex[0], xyzQ);
//	for (int ik = 0 ; ik < kpointIndices.size() ; ++ik )
//	{
//		kandKPPointIndices[ik].first = kpointIndices[ik];
//		kgrid.get_reducible_to_xyz(kpointIndices[ik], xyz);
//		for (int i = 0 ; i < 3 ; ++i)
//			xyz[i] += xyzQ[i];
//		kandKPPointIndices[ik].second = kgrid.get_xyz_to_reducible_periodic(xyz);
//	}
//
//	int nk = kpoints.size()/3;
//
//	elephon::PhononStructure::ElectronPhononCoupling gkkp;
//	elephon::Auxillary::alignedvector::FV gkkpMod2Data;
//
//	// PLEASE NOTE for manual inspection:
//	// The matrix elements computed by generate_gkkp_and_phonon are not symmetrized.
//	// Since the sum of the absolute value squared over the star of a given k point
//	// has the symmetry of the small group of q = k - k', they should be symmetrized
//	// by this small group. In the integral it does not matter, though.
//	gkkp.generate_gkkp_mod_2_of_q(
//			qpointIndex[0],
//			kandKPPointIndices,
//			bandspList, bandsList,
//			resHandl->get_phonon_obj(),
//			resHandl->get_displacement_potential_obj(),
//			resHandl->get_wfct_obj(),
//			gkkpMod2Data);
//	BOOST_REQUIRE_EQUAL(gkkpMod2Data.size(), nk*numBands*numBands*numModes);
//
//	for ( int ik = 0 ; ik < nk ; ++ik )
//	{
//		auto kvect = kgrid.get_vector_direct(kpointIndices[ik]);
//		std::cout << kvect[0] << '\t' << kvect[1] << '\t'<<kvect[2];
//		double integral = 0.0, integralRef = 0.0;
//		for ( int inu = 0 ; inu < numModes ; ++inu)
//		{
//			for ( int ib = 0 ; ib < numBands ; ++ib)
//				for ( int ibp = 0 ; ibp < numBands ; ++ibp)
//					{
//						int cns = ((ik*numBands+ib)*numBands+ibp)*numModes + inu;
//						integral += gkkpMod2Data[cns];
//						integralRef += gkkpDataRef[cns];
//					}
//		}
//		std::cout<<'\t' << integral<<'\t'<<integralRef;
//		std::cout <<'\n';
//	}
//
//	double integral = 0.0, integralRef = 0.0;
//	for ( int ik = 0 ; ik < nk ; ++ik )
//		for ( int ib = 0 ; ib < numBands ; ++ib)
//			for ( int ibp = 0 ; ibp < numBands ; ++ibp)
//				for ( int inu = 0 ; inu < numModes ; ++inu)
//				{
//					int cns = ((ik*numBands+ib)*numBands+ibp)*numModes + inu;
//					integral += gkkpMod2Data[cns];
//					integralRef += gkkpDataRef[cns] ;
//				}
//
//	std::cout << "Integrated |gkkp|^2 for this q = "<< integral << "\t(this code)"
//			<<integralRef <<" (modified espresso)"<<std::endl;
//	BOOST_CHECK_CLOSE(integral, integralRef, 50);
}

BOOST_AUTO_TEST_SUITE_END()
