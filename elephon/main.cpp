/*	This file main.cpp is part of elephon.
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
 *  Created on: 2017
 *      Author: A. Linscheid
 */
#include "IOMethods/Input.h"
#include "ElectronicStructure/LocalDensityOfStates.h"
#include "IOMethods/ChooseInterface.h"

using namespace elephon;

int main(int argc, char* argv[])
{
	IOMethods::Input input(argc,argv);
	auto options = input.get_opts();

	// Open the directory where we expect the preliminary calculation
	auto dataLoader = choose_interface(options);

	// Print the LDOS to file if desired ...
	if ( not options.get_f_ldos().empty() )
	{
		ElectronicStructure::LocalDensityOfStates ldos;
		ldos.compute_ldos(options.get_eldos(), dataLoader);
		ldos.write_file(options.get_f_ldos());
	}

    return 0;
}
