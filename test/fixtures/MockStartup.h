/*	This file MockStartup.h is part of elephon.
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
 *  Created on: Jun 21, 2017
 *      Author: A. Linscheid
 */

#ifndef TEST_FIXTURES_MOCKSTARTUP_H_
#define TEST_FIXTURES_MOCKSTARTUP_H_

#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "IOMethods/InputOptions.h"
#include <string>

namespace elephon
{
namespace test
{
namespace fixtures
{

class MockStartup
{
public:

	boost::filesystem::path get_data_for_testing_dir() const;

	void simulate_elephon_input(
			std::string const & inputFileName,
			std::string const & inputFileContent,
			elephon::IOMethods::InputOptions & inputOpts);

	void write_kpath_file_fcc(std::string filename) const;

	void write_kpath_file_tetra(std::string filename) const;
};

} /* namespace fixtures */
} /* namespace test */
} /* namespace elephon */

#endif /* TEST_FIXTURES_MOCKSTARTUP_H_ */
