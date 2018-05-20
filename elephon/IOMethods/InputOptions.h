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
			lep,
			"Determine if an electron-phonon calculation will be done.",
			"false",
			false,
			bool);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			lp,
			"Determine if a phonon calculation will be done.",
			"false",
			false,
			bool);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			elphd,
			"Where we build the directory structure.\n"
			"Any relative path is relative to root_dir. Only significant if 'lep'  or 'lp' are true",
			"./el_ph/",
			"./el_ph/",
			std::string);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			eld,
			"Where to read/set-up the 'dense' electron states and band structure from/to.\n"
			"If its a relative path, is relative to 'elphd'\n"
			"If empty, the band data from root_dir will used at its place,\n"
			"i.e. no separate 'dense' electronic calculation will be done",
			"electrons",
			"electrons",
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
			"If symDispl is 'true' the code will generate symmetric displacements, i.e. displace in (1,0,0) AND (-1,0,0).\n"
			"This allows to use the second order accurate derivative and leads to more accurate results.\n"
			"If false, the code will only generate (1,0,0) since, mathematically, the derivative is symmetric.",
			"true",
			true,
			bool);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			kdense,
			"Monkhorst-Pack point grid sampling the Brillouin zone.\n"
			"Can be a list of 3 non-negative integers which will be the grid devision.\n"
			"Setting a dimension to 0 means no interpolation in this direction.\n"
			"Can also be single positive integer, in which case it will be a scaling to the read-in grid.\n",
			"4",
			{ 4 },
			std::vector<int>);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			ksdense,
			"Monkhorst-Pack point grid shift in the Brillouin zone.\n"
			"Must be a list of 3 non-negative floating point values smaller 1.0 which\n"
			"will be the shift grid relative to the grid spacing.\n"
			"E.g. '0.5 0.0 0.0' will shift the grid by 0.5/kdense[0] in x direction if \n"
			"kdense determines the grid directly.\n",
			"( 0.0 0.0 0.0 )",
			{ 0.0 COMMA_SUBSTITUTION 0.0 COMMA_SUBSTITUTION 0.0 },
			std::vector<double>);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			kscell,
			"Monkhorst-Pack point grid sampling the Brillouin zone for the\n"
			"supercell calculation. By 'default' it will be the original grid divided by the supercell dimension.\n"
			"Can be a list of 3 non-negative integers which will be the grid devision.\n"
			"Setting a dimension to 0 means 'default' in this direction.\n",
			"( 0 , 0 , 0 )",
			{0 COMMA_SUBSTITUTION 0 COMMA_SUBSTITUTION 0},
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
			ea2f,
			"The energies relative to the Fermi level at which the a2F will be calculated.\n"
			"If more than one number is given, the 'f_a2f' will be preceded by an integer number indicating the position\n"
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
			"Only symmetry non-equivalent extrema are recognized, i.e. the code searches the irreducible zone.\n"
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

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			symOut,
			"Switch to symmetrize output data.\n"
			"Setting this to true means data will be symmetrized using the symmetries of the crystal",
			"true",
			true,
			bool);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			f_a2F,
			"filename of the a2F file.\n"
			"If set, the code will write the Eliashberg function to this location.\n"
			" Default is relative to the 'elphd' option.\n"
			"I.e. the default is '<elphd>/a2F.dat'\n"
			"If empty, the file will not be written.",
			"a2F.dat",
			"a2F.dat",
			std::string);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			f_phdos,
			"filename of the phonon DOS file.\n"
			"If set, the code will compute and write the phonon dos to this location.\n"
			"Default is relative to the 'elphd' option.\n"
			"I.e. the default is '<elphd>/phDOS.dat'\n"
			"If empty, the file will not be written.",
			"phDOS.dat",
			"phDOS.dat",
			std::string);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			phrange,
			"Energy range in THz that the phonons and the a2F will be computed in.\n"
			"Input mode 1 is a single number 'Omega', meaning a range [0,Omega]\n"
			"Input mode 2 are two numbers 'Omega1, Omega2', meaning a range [Omega1,Omega2]\n"
			"Input mode 3 is empty, meaning automatic, i.e. the code will attempt to create\n"
			"\t a range that covers the full spectrum\n",
			"( empty <automatic> )",
			{ },
			std::vector<double>);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			phnpts,
			"Number of sampling points for the phonon DOS and the a2F.\n"
			"Note: For the a2F integration, higher numbers require the Fermi surface sampling to be increased or the data\n"
			"can become noisy.\n",
			"100",
			100,
			int);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			f_dos,
			"If this string is set, the code will compute and print the electronic density of states to this file.\n"
			"The reciprocal grid that will be used to compute it is determined by fftd and ffts.\n"
			"The range will be determined by ewinbnd with a number of sampling points edosnpts.\n"
			"The method is 'tetrahedra'",
			"",
			"",
			std::string);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			edosnpts,
			"Number of sampling points for the electronic DOS.\n",
			"100",
			100,
			int);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			f_kpath,
			"Where to read a band structure path from, both phonons and electronic band structre.\n"
			"Any relative path is relative to root_dir. Only significant if 'f_bands' or 'f_ph_bands' is set\n."
			"\nFormat: <label> kx ky kz <NumPointsSegment> <labelEnd> kEndx kEndy kEndz\n"
			"Labels can be teX code, which will be used in any gnuplot scripts generated by the code.\n"
			"\nExample:\n"
			"$\\Gamma$ 0.0000  0.0000  0.0000 40 X        0.5000  0.0000  0.0000\n"
			"X        0.5000  0.0000  0.0000 40 M        0.5000 -0.5000  0.0000\n",
			"",
			"",
			std::string);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			f_bands,
			"If this string is set, the code will compute and print the electronic band structure to\n"
			"this file. The code will use the fft interpolated mesh starting from the dense band structure\n"
			"and for there use a linear interpolation to obtain the energy eigenvalues at the k points\n"
			"generated according to the file 'f_kpath'\n"
			"Any relative path is relative to root_dir. 'f_kpath' must be set to a file with the k path.\n"
			"The energy window is according to 'ewinbnd'",
			"",
			"",
			std::string);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			f_ph_bands,
			"If this string is set, the code will compute and print the phononic band structure to\n"
			"this file. The code will use the matrix of force constants, fourier transform and compute\n"
			"eigenvalues at the q points generated according to the file 'f_kpath'\n"
			"Any relative path is relative to 'elphd'. 'f_kpath' must be set to a file with the q path.\n"
			"The frequency window is according to 'phrange'",
			"",
			"",
			std::string);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			dvscfc,
			"Defines the coarse gaining of q space on which the dvscf will be caluclated\n"
			"Input mode 1 a single arbitrary negative number in which case the normal k grid will be used\n"
			"Input mode 2 is empty or a single number that is zero where coarse graining is off.\n"
			"Input mode 3 is the explicit grid with 3 integer numbers > 0.\n",
			"( -1 <normal k grid> )",
			{ -1 },
			std::vector<int>);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			EliT,
			"Mode for the isotropic Elishberg equations.\n"
			"The following to options are available:\n"
			"\t'' empty, means do not perform an Elishberg calculation.\n"
			"\t'findTc' where the problem determines the temperature in K where the gap vanishes\n"
			"\t<some list of real numbers> for example '0.1 0.2 0.6 1.5' where the program runs exactly the temperatures specified.\n"
			"The default is 'findTc' if the option f_a2F or Eli_f_a2F has a file path set , and '' otherwise, i.e. if you \n"
			"are running an electron-phonon coupling calculation, elephon will automatically find the Eliashberg Tc for you unless told otherwise.\n",
			"",
			"",
			std::string);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			EliTcThr,
			"Accuracy for the Tc calculation if Tc is searched.\n"
			"The code predicts the Tc from a previous set of runs. It then compares the accuracy\n"
			"of this prediction as it is approaching the critical temperature. \n"
			"If the next prediction is within this window, the code returns the most recent value for Tc\n",
			"0.01",
			0.01,
			double);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			matsC,
			"Cutoff in eV for the Matsubara frequency in an Eliashberg calculation\n"
			"The default is chosen to be conservative for moderate Tc Superconductors.",
			"2.0",
			2.0,
			double);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			muStar,
			"Renormalized isotropic Coulomb interaction at the Fermi level.\n"
			"The default is chosen to be realistic for bulk Superconductors.\n"
			"For multiband systems, this must be of size 1 (taken to be the same among every band)\n"
			"or be of the size of the band matrix. The elements are taken in the order\n"
			"(band 1, band1) (band 1, band2) ... (band 2, band1) ... (band N, bandN)",
			"0.12",
			{ 0.12 },
			std::vector<double>);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			muStarC,
			"Matsubara frequency cutoff in eV for the renormalization term\n"
			"A sensible value is 10 times the boson coupling energy, i.e. the Debye frequency.\n"
			"We except no or one number. If not specified, the program takes 10 times the phonon energy range.",
			"( empty <automatic> )",
			{ },
			std::vector<double>);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			EliCTh,
			"Convergence threshold for the Eliashberg equations.\n"
			"Specify the percentage up to which the gap is allowed to differ when it reproduces itself through the Eliashberg equations\n"
			"before we consider it converged.\n",
			"1e-5",
			1e-5,
			double);


	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			EliCImp,
			"Convergence threshold starting an heuristic convergence acceleration strategy.\n"
			"It appears that especially close to the critical temperature, the shape of the gap function\n"
			"converges much more quickly than its magnitude. We use that to scale by a power of the ratio, to mimic\n"
			"several iterations. The value specifies the sqrt of the square deviations from the previous iteration relative to the gap value. \n",
			"1e-2",
			1e-2,
			double);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			Eli_f_a2F,
			"The file with the isotropic electron-phonon coupling Eliashberg function a2F(w).\n"
			"Format is for each data point frequency w_i [in THz] and then a2F(w_i)\n"
			"the default is f_a2F, if that was set. If it is empty, this value must be specified.\n",
			"( value of f_a2F )",
			"",
			std::string);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			EliMix,
			"Mixing parameter for the Eliashberg equations.\n"
			"Range is from 0 to 1 (excluded) where 0 means an iteration is not averaged with the previous one.\n",
			"0.3",
			0.3,
			double);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			EliThrZero,
			"Threshold what number is considered zero.\n"
			"If the absolute value of any number is below this value, it will not impact convergence\n"
			"If all values of a function, in particular the gap, it is considered zero. Units depend on the context.\n"
			"Thus, the system is considered to be a non-superconductor if every value of the gap drops below this number.\n"
			"For this comparison with the gap, we measure in eV.\n",
			"1e-6",
			1e-6,
			double);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			EliNMax,
			"Maximal number of iterations before the iteration of the Eliashberg equations is aborted.\n"
			"NOTE: Exactly at Tc, the Eliashberg equations become 'margial', i.e. the gap changes at an infinitesimaly small rate.\n"
			"Close to Tc, this number may have to be increased.",
			"10000",
			10000,
			int);
};

} /* namespace IOMethods */
} /* namespace elephon */
#endif /* ELEPHON_IOMETHODS_INPUTOPTIONS_H_ */
