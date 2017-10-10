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
