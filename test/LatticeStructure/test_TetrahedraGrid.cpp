/*	This file test_TetrahedraGrid.cpp is part of elephon.
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
 *  Created on: Oct 8, 2017
 *      Author: A. Linscheid
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE TetraGrids
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "LatticeStructure/TetrahedraGrid.h"
#include "LatticeStructure/RegularSymmetricGrid.h"
#include "fixtures/DataLoader.h"
#include "fixtures/MockStartup.h"
#include <memory>

void output_wire_frame(elephon::LatticeStructure::TetrahedraGrid const & tetra)
{
	int c = 0;
	for ( auto t : tetra.get_tetra_list() )
	{
		std::vector<double> p0, v123;
		t.compute_corner_vectors(p0, v123);

		std::cout << p0[0]<< '\t'<< p0[1]<< '\t' << p0[2]<< '\t'<<v123[0] << '\t' <<v123[1] << '\t' <<v123[2] << '\t' <<c<<'\n';
		std::cout << p0[0]<< '\t'<< p0[1]<< '\t' << p0[2]<< '\t'<<v123[3+0] << '\t' <<v123[3+1] << '\t' <<v123[3+2] << '\t' <<c <<'\n';
		std::cout << p0[0]<< '\t'<< p0[1]<< '\t' << p0[2]<< '\t'<<v123[6+0] << '\t' <<v123[6+1] << '\t' <<v123[6+2] << '\t' <<c<<'\n';

		std::cout << p0[0]+v123[0]<< '\t'<< p0[1]+v123[1]<< '\t' << p0[2]+v123[2]<<
				'\t'<<v123[3+0]-v123[0] << '\t' <<v123[3+1]-v123[1] << '\t' <<v123[3+2]-v123[2] << '\t' <<c<<'\n';
		std::cout << p0[0]+v123[0]<< '\t'<< p0[1]+v123[1]<< '\t' << p0[2]+v123[2]<<
				'\t'<<v123[6+0]-v123[0] << '\t' <<v123[6+1]-v123[1] << '\t' <<v123[6+2]-v123[2] << '\t' <<c<<'\n';
		std::cout << p0[0]+v123[3+0]<< '\t'<< p0[1]+v123[3+1]<< '\t' << p0[2]+v123[3+2] <<
				'\t'<<v123[6+0]-v123[3+0] << '\t' <<v123[6+1]-v123[3+1] << '\t' <<v123[6+2]-v123[3+2]<< '\t' <<c <<'\n';
		c++;
	}
}

BOOST_AUTO_TEST_CASE( Tetrahedra_partitioning_with_sym )
{
	elephon::test::fixtures::DataLoader dl;
	elephon::LatticeStructure::Symmetry sym = dl.create_partial_sym();
	elephon::LatticeStructure::RegularSymmetricGrid g;
	g.initialize({5,5,5}, 1e-6, {0,0,0}, sym);

	elephon::LatticeStructure::TetrahedraGrid tetra;
	tetra.initialize(std::make_shared<decltype(g)>(g));

	int nTetraTotal = 0;
	for ( auto t : tetra.get_tetra_list() )
		nTetraTotal += t.get_multiplicity();
	BOOST_CHECK_EQUAL(g.get_grid_dim()[0]*g.get_grid_dim()[1]*g.get_grid_dim()[2]*6, nTetraTotal);
}

BOOST_AUTO_TEST_CASE( Tetrahedra_partitioning_no_sym )
{
	elephon::LatticeStructure::RegularSymmetricGrid g;
	g.initialize({1,1,1}, 1e-6, {0,0,0});

	elephon::LatticeStructure::TetrahedraGrid tetra;
	tetra.initialize(std::make_shared<decltype(g)>(g));

	BOOST_CHECK_EQUAL( tetra.get_n_tetra() , g.get_grid_dim()[0]*g.get_grid_dim()[1]*g.get_grid_dim()[2]*6 );

	// cross check volume and stuff
	auto vol = [&g] (elephon::LatticeStructure::TetrahedraGrid::Tetrahedra const & t) {
		std::vector<double> p0;
		std::vector<double> v123;
		t.compute_corner_vectors(p0, v123);
		double A11 = v123[3*0+0];double A12 = v123[3*0+1];double A13 = v123[3*0+2];
		double A21 = v123[3*1+0];double A22 = v123[3*1+1];double A23 = v123[3*1+2];
		double A31 = v123[3*2+0];double A32 = v123[3*2+1];double A33 = v123[3*2+2];
		return std::abs(A11*(A22*A33-A23*A32)-A12*(A21*A33-A23*A31)+A13*(A21*A32-A22*A31))/6.0;
	};

	double v0 = 1.0 / static_cast<double>(tetra.get_n_tetra());
	double volSum = 0;
	for ( auto t : tetra.get_tetra_list() )
	{
		double v = vol(t);
		BOOST_CHECK_CLOSE(v0, v, 1e-6 );
		volSum += v;
	}
	BOOST_CHECK_CLOSE(volSum, 1.0, 1e-6);

	// check slightly more complicated grid such that the other diagonal will be chosen for splitting
	g.initialize({2,3,4}, 1e-6, {0,0,0});
	tetra.initialize(std::make_shared<decltype(g)>(g));
	BOOST_CHECK_EQUAL( tetra.get_n_tetra() , g.get_grid_dim()[0]*g.get_grid_dim()[1]*g.get_grid_dim()[2]*6 );
	v0 = 1.0 / static_cast<double>(tetra.get_n_tetra());
	volSum = 0;
	for ( auto t : tetra.get_tetra_list() )
	{
		double v = vol(t);
		BOOST_CHECK_CLOSE(v0, v, 1e-6 );
		volSum += v;
	}
	BOOST_CHECK_CLOSE(volSum, 1.0, 1e-6);
}

BOOST_AUTO_TEST_CASE( Tetrahedra_data_grid )
{
	elephon::LatticeStructure::RegularSymmetricGrid g;
	g.initialize({2,2,2}, 1e-6, {0,0,0});
	elephon::LatticeStructure::TetrahedraGrid tetra;
	tetra.initialize(std::make_shared<decltype(g)>(g));

	for ( auto t : tetra.get_tetra_list() )
	{
		auto ind = t.get_corner_indices();
		std::vector<double> p0123;
		t.compute_corner_points(p0123);
		for ( int i = 0 ; i < 4 ; ++i )
		{
			std::vector<int> xyz(3);
			for ( int j = 0 ; j < 3 ; ++j)
				xyz[j] = std::floor(p0123[i*3+j]*g.get_grid_dim()[j]+0.5);
			auto red = g.get_xyz_to_reducible_periodic(xyz);
			BOOST_CHECK_EQUAL(ind[i], red);
		}
	}

}

BOOST_AUTO_TEST_CASE( Tetrahedra_data_grid_reducible )
{
	elephon::LatticeStructure::RegularSymmetricGrid g;
	elephon::LatticeStructure::RegularSymmetricGrid g_nosym;
	elephon::test::fixtures::DataLoader dl;
	elephon::LatticeStructure::Symmetry sym = dl.create_partial_sym();
	sym.set_reciprocal_space_sym();
	elephon::LatticeStructure::Symmetry id;
	id.set_reciprocal_space_sym();
	g.initialize({1,1,2}, 1e-6, {0,0,0}, sym);
	g_nosym.initialize(g.get_grid_dim(), 1e-6, g.get_grid_shift(), id);
	elephon::LatticeStructure::TetrahedraGrid tetra;
	elephon::LatticeStructure::TetrahedraGrid tetraNoSym;
	tetra.initialize(std::make_shared<elephon::LatticeStructure::RegularSymmetricGrid>(g));
	tetraNoSym.initialize(std::make_shared<elephon::LatticeStructure::RegularSymmetricGrid>(g_nosym));

	// create symmetric data in the irreducible zone
	std::vector<double> dataIr(g.get_np_irred());
	for ( int ikir = 0 ; ikir < g.get_np_irred() ; ++ikir )
	{
		dataIr[ikir] = ikir;
	}
	elephon::LatticeStructure::DataRegularGrid<double> dgrid;
	dgrid.initialize(1, 0.0, dataIr, g);

	// expand to the reducible zone
	elephon::LatticeStructure::DataRegularGrid<double> dgridRed;
	std::vector<double> dataRed;
	dgrid.generate_reducible_data({0}, dataRed);

	BOOST_REQUIRE(dataRed.size() == g.get_np_red());
	dgridRed.initialize(1, 0.0, dataRed, g_nosym);

	// now test if the data indices from the reducible grid
	// match the irreducible one
	for ( auto t : tetraNoSym.get_tetra_list() )
	{
		auto ind = t.get_corner_indices();
		for ( auto & i  : ind )
		{
			i = g.get_maps_red_to_irreducible()[i];
			BOOST_CHECK_EQUAL(double(i), dataIr[i]);
		}
	}

	for ( int itir = 0 ; itir < tetra.get_tetra_list().size(); ++itir )
	{
		auto t = tetra.get_tetra_list()[itir];
		auto star = tetra.get_irreducible_to_reducible(itir);
		for ( auto s : star )
		{
			auto tr = tetra.get_reducible_tetra_list()[s];
			std::vector<double> cornerEnergiesRed(4), cornerEnergiesIrred(4);
			for ( int ic = 0 ;ic < tr.get_corner_indices().size() ; ++ic )
			{
				cornerEnergiesRed[ic] = dgridRed.read(tr.get_corner_indices()[ic], 0);
				cornerEnergiesIrred[ic] = dgrid.read(t.get_corner_indices()[ic], 0);
			}
			std::sort(cornerEnergiesRed.begin(), cornerEnergiesRed.end());
			std::sort(cornerEnergiesIrred.begin(), cornerEnergiesIrred.end());

			for ( int ic = 0 ;ic < tr.get_corner_indices().size() ; ++ic )
				if ( fabs(cornerEnergiesRed[ic] - cornerEnergiesIrred[ic]) > 1e-6)
					std::cout << t.get_corner_indices()[0] << "\t" << t.get_corner_indices()[1] << "\t"
					 << t.get_corner_indices()[2] << "\t" << t.get_corner_indices()[3] << "\t"
					 << tr.get_corner_indices()[0] << "\t" << tr.get_corner_indices()[1] << "\t"
					 << tr.get_corner_indices()[2] << "\t" << tr.get_corner_indices()[3] << "\t"<< std::endl;

			for ( int ic = 0 ;ic < tr.get_corner_indices().size() ; ++ic )
				BOOST_CHECK_CLOSE(cornerEnergiesRed[ic], cornerEnergiesIrred[ic], 1e-6);
		}
	}

	auto mm = dgrid.get_min_max();
	const int numEne = 100;
	std::vector<double> energies(numEne);
	for (int ie = 0 ; ie < numEne; ++ie)
		energies[ie] = mm.first + (mm.second-mm.first)*(ie+0.5)/numEne;


	std::vector<double> irredDos;
	dgrid.compute_DOS_tetra(
			std::make_shared<elephon::LatticeStructure::TetrahedraGrid>(tetra),
			energies,
			irredDos);

	std::vector<double> redDos;
	dgridRed.compute_DOS_tetra(
			std::make_shared<elephon::LatticeStructure::TetrahedraGrid>(tetraNoSym),
			energies,
			redDos);

	BOOST_REQUIRE_EQUAL(redDos.size(), irredDos.size());
}

BOOST_AUTO_TEST_CASE( Tetrahedra_fcc_Al_vasp )
{
	elephon::test::fixtures::MockStartup ms;
	auto rootDir = ms.get_data_for_testing_dir() / "Al" / "vasp" / "fcc_primitive";
	elephon::test::fixtures::DataLoader dl;
	auto resHndler = dl.create_resource_handler(std::string()+
			"root_dir = "+rootDir.string()
			);
	auto g = resHndler->get_electronic_bands_obj()->get_grid();

	elephon::LatticeStructure::TetrahedraGrid tetra;
	tetra.initialize(std::make_shared<decltype(g)>(g));

	double unity = 0;
	int nTetraTotal = 0;
	for ( auto t : tetra.get_tetra_list() )
	{
		unity += double(t.get_multiplicity()) / double(tetra.get_n_reducible_tetra());
		nTetraTotal += t.get_multiplicity();
	}

	BOOST_CHECK_EQUAL(g.get_grid_dim()[0]*g.get_grid_dim()[1]*g.get_grid_dim()[2]*6, nTetraTotal);

	BOOST_CHECK_CLOSE(unity, 1.0, 1e-6);
}

BOOST_AUTO_TEST_CASE( Tetrahedra_fcc_Al_vasp_tetra_dos )
{
	elephon::test::fixtures::MockStartup ms;
	auto rootDir = ms.get_data_for_testing_dir() / "Al" / "vasp" / "fcc_primitive";
	elephon::test::fixtures::DataLoader dl;
	auto resHndler = dl.create_resource_handler(std::string()+
			"root_dir = "+rootDir.string()
			);
	auto bands = resHndler->get_electronic_bands_obj();

	std::vector<int> bndIndices(bands->get_nBnd());
	for (int ibnd = 0 ; ibnd < bndIndices.size(); ++ibnd)
		bndIndices[ibnd] = ibnd;

	std::vector<double> reducibleBndData;
	bands->generate_reducible_data(bndIndices, reducibleBndData);
	auto reducibleBands = std::make_shared<elephon::ElectronicStructure::ElectronicBands>();
	auto gridNoSym = std::make_shared<elephon::LatticeStructure::RegularSymmetricGrid>();
	elephon::LatticeStructure::Symmetry idSym;
	idSym.set_reciprocal_space_sym();
	gridNoSym->initialize( 	bands->get_grid().get_grid_dim(),
							bands->get_grid().get_grid_prec(),
							bands->get_grid().get_grid_shift(),
							idSym,
							bands->get_grid().get_lattice()	);
	reducibleBands->initialize(bndIndices.size(), 0, reducibleBndData, *gridNoSym);

	auto mm = bands->get_min_max();
	const int numEne = 100;
	std::vector<double> energies(numEne);
	for (int ie = 0 ; ie < numEne; ++ie)
		energies[ie] = mm.first + (mm.second-mm.first)*(ie+0.5)/numEne;


	auto tetraReducible = std::make_shared<elephon::LatticeStructure::TetrahedraGrid>();
	tetraReducible->initialize(gridNoSym);

	auto tetra = std::make_shared<elephon::LatticeStructure::TetrahedraGrid>();
	tetra->initialize( std::make_shared<elephon::LatticeStructure::RegularSymmetricGrid>( bands->get_grid() ));

	std::vector<double> irredDos;
	bands->compute_DOS_tetra(tetra, energies, irredDos);

	std::vector<double> redDos;
	reducibleBands->compute_DOS_tetra(tetraReducible, energies, redDos);

	BOOST_REQUIRE_EQUAL(redDos.size(), irredDos.size());

	double numElecRed = 0, numElecIrrRed =0;
	double dE = (mm.second-mm.first)/numEne;
	for ( int ie = 0 ; ie < numEne ;++ie)
	{
		numElecRed += redDos[ie]*dE;
		numElecIrrRed += irredDos[ie]*dE;
	}
	BOOST_CHECK_CLOSE(numElecRed, 5, 1);
	BOOST_CHECK_CLOSE(numElecIrrRed, 5, 1);

	for ( int ie = 0 ; ie < numEne ;++ie)
	{
		BOOST_CHECK_CLOSE(redDos[ie], irredDos[ie], 1e-2);
	}
}
