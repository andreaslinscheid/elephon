/*	This file InputOptions.h is part of elephon.
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
 *  Created on: May 17, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_IOMETHODS_INPUTOPTIONS_H_
#define ELEPHON_IOMETHODS_INPUTOPTIONS_H_

#include "IOMethods/InputBase.h"
#include <vector>
#include <string>

namespace elephon
{
namespace IOMethods
{

class InputOptions : public InputBase<InputOptions>
{
	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			root_dir,
			"Where the initial run took place and we find all initial data.",
			"./",
			"./",
			std::string);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			elphd,
			"Where we build the directory structure.",
			"./el_ph/",
			"./el_ph/",
			std::string);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			scell,
			"Supercell multiplyer in each direction x,y and z. Default is the trivial supercell.",
			"( 1 , 1 , 1 )",
			{1 COMMA_SUBSTITUTION 1 COMMA_SUBSTITUTION 1},
			std::vector<int>);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			numFS,
			"Total number points k points sampling the Fermi surface.",
			"1000",
			1000,
			int);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			magdispl,
			"The magnitude in Angstrom of a displacement in Cartesian coordinates",
			"0.01",
			0.01,
			double);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			symDispl,
			"If symDispl is 'true' the code will generate symmetric displacements, i.e. displace in (1,0,0) AND (-1,0,0)."
			"This allows to use the second order accurate derivative and leads to more accurate results."
			"If false, the code will only generate (1,0,0) since, mathematically, the derivative is symmetric.",
			"true",
			true,
			bool);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			kdense,
			"Monkhorst-Pack point grid sampling the Brillouin zone for the\n"
			"calculation of the electronic structure for the overlap matrix elements.\n"
			"A value <= 0 means automatic which means we scale up by a factor of 5 in each direction",
			"( -1 , -1 , -1 )",
			{-1 COMMA_SUBSTITUTION -1 COMMA_SUBSTITUTION -1},
			std::vector<int>);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			kscell,
			"Monkhorst-Pack point grid sampling the Brillouin zone for the\n"
			"supercell calculation. A value <= 0 means automatic which \n"
			"means we scale down by a factor of 1/scell[i] in each direction",
			"( -1 , -1 , -1 )",
			{-1 COMMA_SUBSTITUTION -1 COMMA_SUBSTITUTION -1},
			std::vector<int>);
};

} /* namespace IOMethods */
} /* namespace elephon */
#endif /* ELEPHON_IOMETHODS_INPUTOPTIONS_H_ */
