/*	This file MockStartup.cpp is part of elephon.
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

#include "fixtures/MockStartup.h"
#include "IOMethods/Input.h"
#include "AtomicSite/AtomSiteData.h"
#include "AtomicSite/SphericalHarmonicExpansion.h"
#include "symmetry/SymmetryOperation.h"
#include "fixtures/DataLoader.h"
#include <fstream>
#include <stdexcept>

namespace elephon
{
namespace test
{
namespace fixtures
{

boost::filesystem::path
MockStartup::get_data_for_testing_dir() const
{
	return boost::filesystem::path(__FILE__).parent_path()  / ".." / "data_for_testing";
}

void
MockStartup::simulate_elephon_input(
		std::string const & inputFileName,
		std::string const & inputFileContent,
		elephon::IOMethods::InputOptions & inputOpts)
{
	//here we create the test input file
	std::ofstream file( inputFileName.c_str() );
	if ( not file.good() )
		throw std::runtime_error(std::string("Can't open file ")+inputFileName+" for opening");
	file << inputFileContent;
	file.close();

	//here we read the input file via elephon's input mechanism
	std::string progn("program name");
	char * prog = new char [progn.size()+1];
	std::copy(progn.begin(),progn.end(),prog);
	prog[progn.size()] = '\0';
	char * arg =  new char [inputFileName.size()+1];
	std::copy(inputFileName.begin(),inputFileName.end(),arg);
	arg[inputFileName.size()] = '\0';
	char *argv[] = {prog, arg, NULL};
	int argc = sizeof(argv) / sizeof(char*) - 1;
	elephon::IOMethods::Input input(argc,argv);
	inputOpts = input.get_opts();
	delete [] prog;
	delete [] arg;

	boost::filesystem::remove( boost::filesystem::path(inputFileName) );
}

void
MockStartup::write_kpath_file_tetra(std::string filename) const
{
	std::string content = ""
	"$\\Gamma$ 0.0000  0.0000  0.0000 40 X        0.5000  0.0000  0.0000\n"
	"X        0.5000  0.0000  0.0000 40 M        0.5000 -0.5000  0.0000\n"
	"M        0.5000 -0.5000  0.0000 57 $\\Gamma$ 0.0000  0.0000  0.0000\n";

	std::ofstream file(filename.c_str());
	if ( ! file.good() )
		throw std::runtime_error("Problem opening file");

	file << content << std::endl;
	file.close();
}

void
MockStartup::write_kpath_file_fcc(std::string filename) const
{
	std::string content = ""
	"$\\Gamma$ 0.0000  0.0000  0.0000 40 X        0.0000  0.5000  0.5000\n"
	"X        0.0000  0.5000  0.5000 40 W        0.2500  0.7500  0.5000\n"
	"W        0.2500  0.7500  0.5000 40 K        0.3750  0.7500  0.3750\n"
	"K        0.3750  0.7500  0.3750 57 $\\Gamma$ 0.0000  0.0000  0.0000\n";

	std::ofstream file(filename.c_str());
	if ( ! file.good() )
		throw std::runtime_error("Problem opening file");

	file << content << std::endl;
	file.close();
}

std::shared_ptr<const AtomicSite::AtomSiteData>
MockStartup::get_mock_AtomSiteData()
{
	if (! this->mockAtomSiteConstantPlusCosX_)
	{
		std::vector<double> pos{0.25, 0.125, 0.0625};
		LatticeStructure::Atom testAtom(20.0, "Bs", pos, {false, false, false});

		AtomicSite::RadialGrid rgrid;
		const int numRadPts = 50;
		std::vector<double> gridPoints(numRadPts);
		const double radius = 0.5;
		for (int ir = 0 ; ir < gridPoints.size(); ++ir)
			gridPoints[ir] = radius*std::pow((0.5+ir)/static_cast<double>(gridPoints.size()), 2);
		rgrid.initialize(testAtom.get_position(), radius, std::move(gridPoints));

		const int lmax = 5;
		Auxillary::alignedvector::ZV angularData(std::pow(lmax+1,2)*numRadPts, std::complex<double>(0));
		AtomicSite::SphericalHarmonicExpansion she;
		she.initialize(lmax, std::move(angularData), rgrid);

		// constant
		for (int ir = 0 ; ir <numRadPts; ++ir)
			she(ir, 0, 0) = M_PI;
		// ~cos(x)
		for (int ir = 0 ; ir <numRadPts; ++ir)
			she(ir, 0, 1) = 2.0*M_PI;

		AtomicSite::FrozenCore fc;
		fc.initialize(-1.0, std::vector<double>(rgrid.get_num_R(),0.0), std::move(rgrid));

		std::shared_ptr<AtomicSite::AtomSiteData> tmp = std::make_shared<AtomicSite::AtomSiteData>();
		tmp->initialize(std::move(testAtom), std::move(she), std::move(fc));
		mockAtomSiteConstantPlusCosX_ = tmp;
	}
	return mockAtomSiteConstantPlusCosX_;
}

symmetry::SymmetryOperation
MockStartup::get_90Deg_rot_about_z_trivial_cell()
{
	DataLoader dl;
	if ( ! symmetry_ )
		symmetry_ = std::make_shared<LatticeStructure::Symmetry>(dl.create_partial_sym());
	return symmetry_->get_sym_op(2);
}


} /* namespace fixtures */
} /* namespace test */
} /* namespace elephon */
