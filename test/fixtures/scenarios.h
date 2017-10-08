/*	This file scenarios.h is part of elephon.
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
 *  Created on: Oct 4, 2017
 *      Author: A. Linscheid
 */

#ifndef TEST_FIXTURES_SCENARIOS_H_
#define TEST_FIXTURES_SCENARIOS_H_

#include "IOMethods/ResourceHandler.h"
#include "fixtures/MockStartup.h"
#include "fixtures/DataLoader.h"
#include <boost/filesystem.hpp>
#include <string>
#include <memory>

namespace elephon
{
namespace test
{
namespace fixtures
{
namespace scenarios
{

std::shared_ptr<IOMethods::ResourceHandler>
load_Al_fcc_primitive_vasp_sc2x2x2()
{
	test::fixtures::MockStartup ms;
	auto rootDir = ms.get_data_for_testing_dir() / "Al" / "vasp" / "fcc_primitive" / "phonon_run" ;
	auto phononDir = rootDir / "phonon";
	std::string content = std::string()+
			"scell=2 2 2\n"
			"root_dir="+rootDir.string()+"\n"
			"elphd="+phononDir.string()+"\n"
			"eld="+ (ms.get_data_for_testing_dir() / "Al" / "vasp" / "fcc_primitive").string()+"\n"
			"f_a2F="+ (ms.get_data_for_testing_dir() / "Al" / "vasp" / "fcc_primitive" / "output" / "a2F.dat").string()+"\n"
			"phrange = 10\n"
			"numFS = 10"
			"";
	test::fixtures::DataLoader dl;
	return dl.create_resource_handler( content );
}

} /* namespace scenarios */
} /* namespace fixtures */
} /* namespace test */
} /* namespace elephon */

#endif /* TEST_FIXTURES_SCENARIOS_H_ */
