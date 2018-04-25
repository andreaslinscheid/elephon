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
namespace AtomicSite {class AtomSiteData;};
namespace symmetry { class SymmetryOperation; };
namespace LatticeStructure { class Symmetry; } ;

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

	/**
	 * Obtain sample data of an atom at site (0.25 0.125 0.0625) with
	 * max l quantum number 5. Only the constant and the ~cos(x) channel
	 * is set to PI and 2PI, respectively. The number of radial points is 50
	 * and the radius of that data 0.5
	 *
	 * @return	shared ptr with the data
	 */
	std::shared_ptr<const AtomicSite::AtomSiteData> get_mock_AtomSiteData();

	/**
	 * Return the symmetry operation that rotates by 90 degrees about the z axis
	 * for a trivial unit cell with extend of 1 in each direction.
	 *
	 * @return the symmetry operation.
	 */
	symmetry::SymmetryOperation get_90Deg_rot_about_z_trivial_cell();

private:

	std::shared_ptr<const AtomicSite::AtomSiteData> mockAtomSiteConstantPlusCosX_;

	std::shared_ptr<LatticeStructure::Symmetry> symmetry_;
};

} /* namespace fixtures */
} /* namespace test */
} /* namespace elephon */

#endif /* TEST_FIXTURES_MOCKSTARTUP_H_ */
