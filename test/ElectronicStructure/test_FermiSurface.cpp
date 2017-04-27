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
#include "ElectronicStructure/FermiSurface.h"

BOOST_AUTO_TEST_CASE( Check_DOS_LiFeAs )
{
	//load data from the LiFeAs electronic structure
	std::string LiFeAs_elstr_f = "../ElectronicStructure/LiFeAs_energies.dat";
	std::ifstream input( LiFeAs_elstr_f.c_str(), std::ios::binary );
	std::vector<char> buffer( std::istreambuf_iterator<char>(input),
	            std::istreambuf_iterator<char>());
	double const * data_view = reinterpret_cast<double>(buffer.data());
	size_t nOrb = std::round(data_view[0]);
	size_t nkx = std::round(data_view[1]);
	size_t nky = std::round(data_view[2]);
	size_t nkz = std::round(data_view[3]);
	assert( data_view[0]*data_view[1]*data_view[2])

	elephon::ElectronicStructure::FermiSurface fs;

	BOOST_CHECK( input.get_numFS() == 2000 );
}
