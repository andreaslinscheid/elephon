/*	This file Input.h is part of elephon.
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

#ifndef ELEPHON_IOMETHODS_INPUT_H_
#define ELEPHON_IOMETHODS_INPUT_H_

#include "IOMethods/InputOptions.h"

namespace elephon
{
namespace IOMethods
{

/**
 * The main program configuration is obtained via this class.
 *
 * It reads an input file and allows to give constant read access to the resulting program options.
 */
class Input
{
public:

	/**
	 * Based on program startup options from the command line, obtain the configuration options.
	 *
	 * @param argc	Number of arguments passed via the command line.
	 * @param argv	The arguments passed to the command line.
	 */
	Input( int argc, char* argv[] );

	/**
	 * Obtain read acces to the internally stored program options.
	 *
	 * @return	A constant reference to the InputOptions class.
	 */
	InputOptions const & get_opts() const;

private:

	InputFile inputFile_;

	InputOptions opts_;
};

} /* namespace IOMethods */
} /* namespace elephon */

#endif /* ELEPHON_IOMETHODS_INPUT_H_ */
