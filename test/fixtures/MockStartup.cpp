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
#include <fstream>
#include <stdexcept>

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

} /* namespace fixtures */
} /* namespace test */
