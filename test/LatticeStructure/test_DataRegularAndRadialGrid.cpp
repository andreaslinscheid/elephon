/*	This file test_DataRegularAndRadialGrid.cpp is part of elephon.
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
 *  Created on: Feb 28, 2018
 *      Author: A. Linscheid
 */
#include <boost/test/unit_test.hpp>
#include "LatticeStructure/DataRegularAndRadialGrid.h"
#include "LatticeStructure/Atom.h"
#include "Auxillary/AlignedVector.h"
#include "AtomicSite/AtomSiteData.h"
#include "AtomicSite/SphericalHarmonicExpansion.h"
#include "fixtures/MockStartup.h"
#include "AtomicSite/FrozenCore.h"

BOOST_AUTO_TEST_SUITE( DataRegularAndRadialGrid )

template<typename T> void test_basic();

BOOST_AUTO_TEST_CASE( test_basic_functionality )
{
	// test the basic interface functionality
	test_basic<float>();
	test_basic<std::complex<float>>();
	test_basic<double>();
	test_basic<std::complex<double>>();
}

template<typename T>
void test_basic()
{
	using namespace elephon;
	LatticeStructure::DataRegularAndRadialGrid<T> d;
	elephon::test::fixtures::MockStartup ms;

	T valueRegularGrid = T(1.0);
	Auxillary::alignedvector::aligned_vector<T> regularGridData(125, valueRegularGrid);
	std::vector<AtomicSite::AtomSiteData> radialGridData(1);
	const int atomIndex = 0;
	auto atom = ms.get_mock_AtomSiteData()->get_atom();
	atom.set_position({0,0,0});
	auto dataAtom = ms.get_mock_AtomSiteData()->get_potential_data();
	dataAtom.set_center(atom.get_position());
	radialGridData[atomIndex].initialize(atom, dataAtom, elephon::AtomicSite::FrozenCore());
	d.initialize(LatticeStructure::RegularBareGrid({5,5,5}), regularGridData, radialGridData);

	BOOST_CHECK_EQUAL(d.get_max_num_radial_elements(), 50);
	const int lmax = 5;
	BOOST_CHECK_EQUAL(d.get_max_num_angular_moment_channels(), std::pow(lmax+1,2));
	BOOST_CHECK_EQUAL(d.get_max_angular_moment(), lmax);

	auto d_rot = d;
	d_rot.transform(ms.get_90Deg_rot_about_z_trivial_cell());

	// check that the data we have initialized is where it is supposed to be.
	// NOTE: the data layout is build into the algorithms, thus this must be checked!
	for (auto begin_reg = d.begin_regular_data(); begin_reg != d.end_regular_data(); ++begin_reg)
	{
		auto diff = std::abs(*begin_reg - valueRegularGrid);
		BOOST_CHECK_SMALL(diff, decltype(diff)(1e-8));
	}

	auto begin_rad = d.begin_radial_data(atomIndex);
	for (auto it = begin_rad ; it != d.end_radial_data(atomIndex); ++it )
	{
		const int dist = std::distance(begin_rad, it);
		auto is_in_channel = [&] (int distance, int l, int m) {
			const int lm = Auxillary::memlayout::angular_momentum_layout(l, m);
			return (dist >= lm*d.get_max_num_radial_elements())
					and (dist < (lm+1)*d.get_max_num_radial_elements());
		};
		if ( is_in_channel(dist, /*l=*/0, /*m=*/0) ) // constant range
		{
			BOOST_CHECK_SMALL(std::abs(*it-M_PI), decltype(std::abs(*it))(1e-8));
		}
		else if (is_in_channel(dist, /*l=*/1, /*m=*/0)) // cos(x)
		{
			BOOST_CHECK_SMALL(std::abs(*it-2.0*M_PI), decltype(std::abs(*it))(1e-8));
		}
		else
			BOOST_CHECK_SMALL(std::abs(*it), decltype(std::abs(*it))(1e-8));
	}

	BOOST_CHECK_EQUAL(d.get_num_atom_data_sets(), 1);
};

BOOST_AUTO_TEST_SUITE_END()
