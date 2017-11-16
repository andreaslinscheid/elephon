/*	This file test_AlphaSquaredF.cpp is part of elephon.
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
 *  Created on: Sep 30, 2017
 *      Author: A. Linscheid
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE PhononStructure
#include <boost/test/unit_test.hpp>
#include "PhononStructure/AlphaSquaredF.h"
#include "IOMethods/ResourceHandler.h"
#include "fixtures/MockStartup.h"
#include "fixtures/DataLoader.h"
#include "fixtures/scenarios.h"
#include <vector>

BOOST_AUTO_TEST_CASE( a2F_write_file_vasp_fcc_primitive )
{
	//TODO figure out a way of testing the a2F generation that is acceptably fast.
//	auto resourceHandler = elephon::test::fixtures::scenarios::load_Al_fcc_primitive_vasp_sc2x2x2();
//
//	elephon::PhononStructure::AlphaSquaredF a2F;
//	a2F.compute_a2F(resourceHandler);
//
//	auto a2FFilename = boost::filesystem::path(resourceHandler->get_optns().get_f_a2F());
//	boost::filesystem::remove(a2FFilename);
//	a2F.write_a2F_file( a2FFilename.string() );
//	boost::filesystem::remove(a2FFilename);
}
