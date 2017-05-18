/*	This file Input.cpp is part of elephon.
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
 *  Created on: Apr 24, 2017
 *      Author: A. Linscheid
 */

#include "IOMethods/Input.h"
#include <stdexcept>
#include <string>
#include <fstream>
#include <iostream>

namespace elephon
{
namespace IOMethods
{

Input::Input( int argc, char* argv[] )
{
	if ( argc != 2)
		throw std::runtime_error( "\tPass input file as parameter, exiting ... !" );
	inputFile_.read_input_file(argv[1]);

	//Fetch the general configuration variables from the file
	opts_.parse_variables(inputFile_);

	//Flag a warning for all keys that have not been used, as this indicates unintended input
	std::vector<std::string> const unusedKeys = inputFile_.get_list_unread_input_parameters();;
	for ( auto & key : unusedKeys )
		std::cout << "WARNING: Key '"+ key +"' defined in the input file has not been\n"
				" matched to a variable used by the code. Typo?" ;

	//TODO Here could be some consistency checks
}

InputOptions const & Input::get_opts() const
{
	return opts_;
}

} /* namespace IOMethods */
} /* namespace elephon */
