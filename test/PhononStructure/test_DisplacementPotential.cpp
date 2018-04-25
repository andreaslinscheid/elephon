/*	This file test_DisplacementPotential.cpp is part of elephon.
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
 *  Created on: Jul 1, 2017
 *      Author: A. Linscheid
 */
#include <boost/test/unit_test.hpp>
#include "fixtures/MockStartup.h"
#include "fixtures/DataLoader.h"
#include "PhononStructure/DisplacementPotential.h"
#include "PhononStructure/Phonon.h"
#include "LatticeStructure/PrimitiveToSupercellConnection.h"
#include "LatticeStructure/AtomDisplacementCollection.h"
#include "PhononStructure/PotentialChangeIrredDisplacement.h"
#include "fixtures/scenarios.h"
#include <iomanip>

BOOST_AUTO_TEST_SUITE( DisplacementPotential )

BOOST_AUTO_TEST_CASE( Mimimal_example)
{
	using namespace elephon;
	/*
	 * Here we create a displacement potential for a minimal system
	 * where the result is analytically known. The ingredients to this
	 * system are discussed below:
	 */
	// Symmetry
	test::fixtures::DataLoader dl;
	auto symmetry4FoldWithInversion = dl.create_4fold_symmetry();

	// Atom compatible with the symmetry
	LatticeStructure::Atom atom(20.0, "Bs", {0.0, 0.0, 0.0}, {false, false, false});

	// With a trivial cubic lattice of length one, we initialize the unit cell
	auto pc = std::make_shared<LatticeStructure::UnitCell>();
	pc->initialize(std::vector<LatticeStructure::Atom>(1, atom), LatticeStructure::LatticeModule(), symmetry4FoldWithInversion);

	// A 2x2x2 supercell and the relation between the two
	auto sc = std::make_shared<LatticeStructure::UnitCell>(pc->build_supercell(2, 2, 2));
	auto ptos = std::make_shared<LatticeStructure::PrimitiveToSupercellConnection>();
	ptos->initialize(pc, sc);

	// Inititalize the displacement to given magnitude
	const double displacementMagnitude = 0.001;
	auto displ = std::make_shared<LatticeStructure::AtomDisplacementCollection>();
	displ->initialize(pc, true, displacementMagnitude);

	// the regular grids
	// PLEASE NOTE: For even grids in the primitive cell we
	// are creating unsymmetric data that causes errors on the order
	// of ration of the number of surface points / total number of points
	LatticeStructure::RegularBareGrid pcGrid({11, 11, 11});
	LatticeStructure::RegularBareGrid scGrid({22, 22, 22});

	// Carefully set the ground state potential and the displaced potential
	// so that we arrive at a controllable result.

	// Ground state is first:
	// We use a Gaussian potential centered around r=0 on the regular and
	// a -sin(pi*r)/|r| potential on the radial grid.
	Auxillary::alignedvector::DV regularGridData(pcGrid.get_num_points());
	auto gridVectors = pcGrid.get_all_vectors_grid();
	for (int ig = 0 ; ig < pcGrid.get_num_points(); ++ig)
	{
		double x = gridVectors[ig*3+0];
		double y = gridVectors[ig*3+1];
		double z = gridVectors[ig*3+2];
		regularGridData[ig] = std::exp(-(x*x+y*y+z*z))/2.0;
	}
	std::vector<AtomicSite::AtomSiteData> atomData(1);
	AtomicSite::SphericalHarmonicExpansion she;
	AtomicSite::RadialGrid rgrid;
	const double radius = 0.5;
	std::vector<double> radialPoints(31);
	for (int ir = 0 ; ir < radialPoints.size(); ++ir)
		radialPoints[ir] = 0.0001+(radius-0.0001)/(radialPoints.size()-1)*ir;
	rgrid.initialize(atom.get_position(), radius, std::move(radialPoints));
	Auxillary::alignedvector::ZV radialData(rgrid.get_num_R()*36, std::complex<double>(0));
	she.initialize(5, radialData, rgrid);
	for (int ir = 0 ; ir < rgrid.get_num_R(); ++ir)
		she(ir, 0, 0) = -(std::sin(M_PI*rgrid.get_radius(ir)*2.0)/rgrid.get_radius(ir))
							*(0.25-std::pow(rgrid.get_radius(ir),2))*std::sqrt(4.0*M_PI); // -sin(pi|r|)/|r|*(0.5^2-r^2)
	atomData[0].initialize(atom, she);
	LatticeStructure::DataRegularAndRadialGrid<double> gsData;
	gsData.initialize(pcGrid,regularGridData, atomData);

	// Displaced state is next,
	auto d = displ->get_irreducible_displacements()[0].second[0];
	auto pos = d.get_position();
	ptos->primitive_to_supercell_coordinates(pos);
	d.set_position(pos);
	auto displAtom = sc->get_atoms_list()[ptos->primitive_to_supercell_atom_index(0)];
	displAtom.apply_displacement(d);

	LatticeStructure::DataRegularAndRadialGrid<double> dsData;
	regularGridData = Auxillary::alignedvector::DV(scGrid.get_num_points());
	gridVectors = scGrid.get_all_vectors_grid();
	ptos->supercell_to_primitive_coordinates(gridVectors);
	for (int ig = 0 ; ig < scGrid.get_num_points(); ++ig)
	{
		using Algorithms::helperfunctions::nint;
		double x = gridVectors[ig*3+0];
		double y = gridVectors[ig*3+1];
		double z = gridVectors[ig*3+2];
		if (((x>=-0.5) && (x<0.5)) and ((y>=-0.5) && (y<0.5)) and ((z>=-0.5) && (z<0.5)))
		{
			// displace the gaussian potential only in the first cell
			x -= d.get_direction()[0]*d.get_magnitude();
			y -= d.get_direction()[1]*d.get_magnitude();
			z -= d.get_direction()[2]*d.get_magnitude();
		}
		else
		{
			// map back to the first unit cell
			x -= nint(x);
			y -= nint(y);
			z -= nint(z);
		}
		regularGridData[ig] = std::exp(-(x*x+y*y+z*z))/2.0;
	}
	std::vector<AtomicSite::AtomSiteData> atomDataSC(sc->get_atoms_list().size());
	for (int iaSC = 0 ; iaSC < sc->get_atoms_list().size(); ++iaSC)
	{
		if ( iaSC == ptos->primitive_to_supercell_atom_index(0) )
		{
			she.set_center(displAtom.get_position());
			atomDataSC[iaSC].initialize(displAtom, she);
		}
		else
		{
			she.set_center(sc->get_atoms_list()[iaSC].get_position());
			atomDataSC[iaSC].initialize(sc->get_atoms_list()[iaSC], she);
		}
	}
	dsData.initialize(scGrid,regularGridData, atomDataSC);

	// Finally, set the potential change
	PhononStructure::PotentialChangeIrredDisplacement potChange;
	potChange.initialize(displ->get_irreducible_displacements()[0].second[0],
			gsData, dsData, std::make_shared<LatticeStructure::RegularBareGrid>(pcGrid),
			 std::make_shared<LatticeStructure::RegularBareGrid>(scGrid), ptos);
	std::vector<std::shared_ptr<const PhononStructure::PotentialChangeIrredDisplacement>> potChangeV;
	potChangeV.push_back(std::make_shared<PhononStructure::PotentialChangeIrredDisplacement>(std::move(potChange)));

	PhononStructure::DisplacementPotential displacementPot;
	displacementPot.initialize(pc, sc, displ, ptos, pcGrid, scGrid, potChangeV);

	// Regular Grid:
	// We expect the derivative of the Gaussian potential on the regular grid.
	Auxillary::Multi_array<float,3> const & displDataRegular = displacementPot.get_data_regular_grid();
	double diff = 0.0;
	std::vector<int> xyz(3,0);
	// We expect contributions only for the first zone iR = 0
	gridVectors = pcGrid.get_all_vectors_grid();
	for (int ig = 0 ; ig < pcGrid.get_num_points(); ++ig)
	{
		int iR = 0;
		using Algorithms::helperfunctions::nint;
		double x = gridVectors[ig*3+0];
		double y = gridVectors[ig*3+1];
		double z = gridVectors[ig*3+2];
		// in x direction
		double referenceValue = -x*std::exp(-(x*x+y*y+z*z));
		int mu = 0;
		diff += std::abs(referenceValue-displDataRegular[iR][mu][ig]);
		// in y direction
		referenceValue = -y*std::exp(-(x*x+y*y+z*z));
		mu = 1;
		diff += std::abs(referenceValue-displDataRegular[iR][mu][ig]);
		// in z direction
		referenceValue = -z*std::exp(-(x*x+y*y+z*z));
		mu = 2;
		diff += std::abs(referenceValue-displDataRegular[iR][mu][ig]);
	}
	diff /= (displacementPot.get_num_modes()*pcGrid.get_num_points());
	BOOST_CHECK_SMALL(diff, 1e-6);

	// All but the first zone contributions must be zero
	diff = 0.0;
	for (int iR = 1 ; iR < displacementPot.get_num_R(); ++iR)
		for (int mu = 0 ; mu < displacementPot.get_num_modes(); ++mu)
			for (int ig = 0 ; ig < pcGrid.get_num_points(); ++ig)
			{
				diff += std::abs(displDataRegular[iR][mu][ig]);
			}
	BOOST_CHECK_SMALL(diff, 1e-6);

	// Radial grid:
	// We expect contributions only for the first zone iR = 0 where the function
	// is supposed to be
	//		- pi * x,y,z cos(pi*|r|) / |r|^2 + x,y,z * sin(pi*|r|)/ |r|^3
	// where x,y,z is the unit vector in either x, y or z direction.
	// We check this data on a regular grid about the atomic location
	Auxillary::Multi_array<std::complex<float>,4> const & displDataRadial = displacementPot.get_data_radial_grid();
	const int radialDataSizePerLatticeVector = displacementPot.get_num_modes()*displDataRadial.shape()[2]*displDataRadial.shape()[3];
	diff = 0.0;
	for (int iR = 1 ; iR < displacementPot.get_num_R(); ++iR)
		for (int mu = 0 ; mu < displacementPot.get_num_modes(); ++mu)
			for (int lm = 0 ; lm < displDataRadial.shape()[2]; ++lm)
				for (int irad = 0 ; irad < displDataRadial.shape()[3]; ++irad)
					diff += std::abs(displDataRadial[iR][mu][lm][irad]);
	diff /= (displacementPot.get_num_R()-1)*radialDataSizePerLatticeVector;
	BOOST_CHECK_SMALL(diff, 1e-6);

	// define the grid where we verify the data
	diff = 0.0;
	const int nPtsCheckGrid = 30;
	const double minX = rgrid.get_center()[0]-rgrid.get_range_of_definition();
	const double minY = rgrid.get_center()[1]-rgrid.get_range_of_definition();
	const double minZ = rgrid.get_center()[2]-rgrid.get_range_of_definition();
	auto setXYZ_within_radius = [&] (int ix, int iy, int iz, double &x, double &y, double &z) {
		x = minX + 2.0*(ix*rgrid.get_range_of_definition())/nPtsCheckGrid;
		y = minY + 2.0*(iy*rgrid.get_range_of_definition())/nPtsCheckGrid;
		z = minZ + 2.0*(iz*rgrid.get_range_of_definition())/nPtsCheckGrid;
		double r = std::sqrt(	  std::pow(x-rgrid.get_center()[0],2)
					+ std::pow(y-rgrid.get_center()[1],2)
					+ std::pow(z-rgrid.get_center()[2],2));
		if ( (r > 1e-4) && (r < rgrid.get_range_of_definition()) )
			return true;
		return false;
	};
	std::vector<double> atomGridCheck;
	double x, y, z;
	for (int ix = 0 ; ix < nPtsCheckGrid; ++ix)
		for (int iy = 0 ; iy < nPtsCheckGrid; ++iy)
			for (int iz = 0 ; iz < nPtsCheckGrid; ++iz)
				if ( setXYZ_within_radius(ix,iy,iz,x,y,z) )
					atomGridCheck.insert(atomGridCheck.end(), {x, y, z});

	// we start with the x direction
	Auxillary::alignedvector::ZV displPotExpansionData(&displDataRadial[0][0][0][0], &displDataRadial[0][1][0][0]);
	Auxillary::alignedvector::ZV outputDataX, outputDataY, outputDataZ;
	AtomicSite::SphericalHarmonicExpansion shX;
	shX.initialize(she.get_l_max(), displPotExpansionData, rgrid);
	shX.interpolate(atomGridCheck, outputDataX);
	BOOST_REQUIRE_EQUAL(outputDataX.size(), atomGridCheck.size()/3);

	// next is the y direction
	std::copy(&displDataRadial[0][1][0][0], &displDataRadial[0][2][0][0], displPotExpansionData.data());
	AtomicSite::SphericalHarmonicExpansion shY;
	shY.initialize(she.get_l_max(), displPotExpansionData, rgrid);
	shY.interpolate(atomGridCheck, outputDataY);
	BOOST_REQUIRE_EQUAL(outputDataY.size(), atomGridCheck.size()/3);

	// and finally the z direction
	std::copy(&displDataRadial[0][2][0][0], &displDataRadial[0][2][0][0]+displPotExpansionData.size(), displPotExpansionData.data());
	AtomicSite::SphericalHarmonicExpansion shZ;
	shZ.initialize(she.get_l_max(), displPotExpansionData, rgrid);
	shZ.interpolate(atomGridCheck, outputDataZ);
	BOOST_REQUIRE_EQUAL(outputDataZ.size(), atomGridCheck.size()/3);

	int c = 0;
	for (int ix = 0 ; ix < nPtsCheckGrid; ++ix)
	{
		for (int iy = 0 ; iy < nPtsCheckGrid; ++iy)
		{
			for (int iz = 0 ; iz < nPtsCheckGrid; ++iz)
			{
				setXYZ_within_radius(ix,iy,iz,x,y,z);
				if ( setXYZ_within_radius(ix,iy,iz,x,y,z) )
				{
					auto driv = [] (double x, double y, double z, double d) {
						return (-2*M_PI*d*(0.25 - std::pow(x,2) - std::pow(y,2) - std::pow(z,2))*
							      std::cos(2*M_PI*std::sqrt(std::pow(x,2) + std::pow(y,2) + std::pow(z,2))))/
							    (std::pow(x,2) + std::pow(y,2) + std::pow(z,2)) +
							   (d*(0.25 - std::pow(x,2) - std::pow(y,2) - std::pow(z,2))*
							      std::sin(2*M_PI*std::sqrt(std::pow(x,2) + std::pow(y,2) + std::pow(z,2))))/
							    std::pow(std::pow(x,2) + std::pow(y,2) + std::pow(z,2),1.5) +
							   (2*d*std::sin(2*M_PI*std::sqrt(std::pow(x,2) + std::pow(y,2) + std::pow(z,2))))/
							    std::sqrt(std::pow(x,2) + std::pow(y,2) + std::pow(z,2));
					};
					double referenceValue = driv(x,y,z,x);
					diff += std::abs(std::real(outputDataX[c])-referenceValue);
					diff += std::abs(std::imag(outputDataX[c]));

					referenceValue = driv(x,y,z,y);
					diff += std::abs(std::real(outputDataY[c])-referenceValue);
					diff += std::abs(std::imag(outputDataY[c]));

					referenceValue = driv(x,y,z,z);
					diff += std::abs(std::real(outputDataZ[c])-referenceValue);
					diff += std::abs(std::imag(outputDataZ[c]));
					++c;
				}
			}
		}
	}

	diff /= atomGridCheck.size()/3;
	// For the small number of radial points that we use, this result is not
	// very accurate. Since re-fitting is a rather expensive operation and essentially O(n_radial^2),
	// we don't run a more approprite number of points.
	BOOST_CHECK_SMALL(diff, 0.01);
}
//
//BOOST_AUTO_TEST_CASE( build_Al_fcc_primitive )
//{
//	auto resourceHandler = elephon::test::fixtures::scenarios::load_Al_fcc_primitive_vasp_sc2x2x2();
//	auto dvscf = resourceHandler->get_displacement_potential_obj();
//
//	BOOST_REQUIRE( dvscf->get_num_modes() == 3 );
//
//	BOOST_REQUIRE( dvscf->get_num_R() == 3*3*3 );
//
//	//Write the real space variant
//	auto rootDir = boost::filesystem::path( resourceHandler->get_optns().get_root_dir() );
//	boost::filesystem::path dvscfFile = rootDir / "dvscf.dat";
//	dvscf->write_dvscf(0,0,dvscfFile.string());
//
//	//Write the q displacement variant
//	std::vector<double> qVect{ 0.0,0.0,0.0 , 0.25,0.0,0.0, 0.5,0.0,0.0 };
//	std::vector<int> modes{0,1};
//	elephon::Auxillary::Multi_array<double,2> w;
//	elephon::Auxillary::Multi_array<std::complex<double>,3> dynMat;
//	auto ph = resourceHandler->get_phonon_obj();
//	ph->compute_at_q( qVect, w, dynMat );
//
//	boost::filesystem::path dvscfqFile = rootDir / "dvscf_q.dat";
//	dvscf->write_dvscf_q(qVect, modes, w, dynMat, ph->get_masses(), dvscfqFile.string());
//
//	//At this point we can perform tests on the files and its content.
//
//	//Outcomment the following for manual inspection of the files generated
//	BOOST_REQUIRE( boost::filesystem::is_regular_file( dvscfFile ) );
//	boost::filesystem::remove( dvscfFile );
//
//	for ( auto mu : modes )
//	{
//		for ( int iq = 0 ; iq < qVect.size()/3; ++iq)
//		{
//			std::string filename = std::string("dvscf_q_")+std::to_string(iq)+"_"+std::to_string(mu)+".dat" ;
//			BOOST_REQUIRE( boost::filesystem::is_regular_file( rootDir / filename ) );
//			boost::filesystem::remove( rootDir / filename );
//		}
//	}
//}
BOOST_AUTO_TEST_SUITE_END()
