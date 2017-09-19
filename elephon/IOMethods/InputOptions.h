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

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			gPrec,
			"The grid precision which decides component-wise up to which value two grid vectors are considered equal.\n",
			"1e-6",
			1e-6,
			double);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			f_ldos,
			"If this string is set, the code will compute and print the local density of states to this file.\n"
			"Default format is the code that produced the input data.",
			"",
			"",
			std::string);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			eldos,
			"The energies relative to the Fermi level at which the ldos will be calculated.\n"
			"If more than one number is given, the 'f_ldos' will be preceded by an integer number indicating the position\n"
			"In the list.",
			"0.0",
			{ 0.0 },
			std::vector<double>);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			fftd,
			"Grid dimension for the electronic structure fft interpolation.\n"
			"Can be a list of 3 non-negative integers which will be the grid.\n"
			"Setting a dimension to 0 means no interpolation in this direction.\n"
			"Can also be single positive integer, in which case it will be a scaling to the read-in grid.\n",
			"3",
			{ 3 },
			std::vector<int>);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			ffts,
			"Grid shift for the electronic structure fft interpolation.\n"
			"Must be a list of 3 non-negative floating point values smaller 1.0 which"
			"will be the shift grid relative to the grid spacing.\n"
			"E.g. '0.5 0.0 0.0' will shift the grid by 0.5/ffd[0] in x direction if fftd determines the grid directly.\n",
			"0.0 0.0 0.0",
			{0.0 COMMA_SUBSTITUTION 0.0 COMMA_SUBSTITUTION 0.0},
			std::vector<double>);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			f_mtens,
			"If this string is set, the code will compute and print the mass tensor at band extrema to this file.\n"
			"The file is binary and contains 16 floats per band extrema. Per extrema, the code prints:\n"
			"\t1.) '+/-1.0', is the flag for minima (<0) or maxima (>0)\n"
			"\t2.) '<integer as float>', is the band index starting to count from 0\n"
			"\t3.) '<3 floats>', the k location of the maximum in the first BZ in units of the reciprocal lattice\n"
			"\t4.) '3x<1+3 floats>', 3 combinations of (eigenvalue plus eigenvector;"
			"3 components in units of the reciprocal lattice) in ascending order of the eigenvalue\n\n"
			"Only symmetry non-equivalent extrema are recognized\n"
			"",
			"",
			"",
			std::string);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			ewinbnd,
			"Energy window for the band structure. Only bands will be considered that cross this energy window in "
			"some points of the unit cell.\n"
			"Must be empty or a list of 2 floats representing the lower and upper limit in eV, respectively.",
			" <empty> ",
			{ },
			std::vector<double>);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			dmeth,
			"string that sets the way the code computes derivatives of the band structure.\n"
			"Methods are 'fft' or 'pol'. fft requires the bands to be ordered by continuity and not by energy.\n"
			"Currently, this ordering is not implemented.\n"
			"pol perform a local least square fit to a quadratic polynomial using nearest and next nearest neighbor points.",
			"pol",
			"pol",
			std::string);
};

} /* namespace IOMethods */
} /* namespace elephon */
#endif /* ELEPHON_IOMETHODS_INPUTOPTIONS_H_ */
