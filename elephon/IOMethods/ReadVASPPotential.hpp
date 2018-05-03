/*	This file ReadVASPPotential.hpp is part of elephon.
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
 *  Created on: Apr 26, 2018
 *      Author: A. Linscheid
 */

#include "IOMethods/ReadVASPPotential.h"
#include "Algorithms/helperfunctions.hpp"
#include "Auxillary/memory_layout_functions.hpp"
#include <iostream>
#include <stdexcept>
#include <cstring>

namespace elephon
{
namespace IOMethods
{

template<class VTReg, class VTRad, class VTChg>
void
ReadVASPPotential::read_potential_file(
		std::array<int,3> & regularGridDim,
		VTReg & regularData,
		std::vector<int> & angularLMaxPerAtom,
		std::vector<std::vector<double>> & radialPointsPerAtom,
		std::vector<double> & radiusPerAtom,
		std::vector<std::array<double,3>> & centerAtom,
		std::vector<std::string> & namesAtom,
		std::vector<double> & massAtom,
		std::vector<VTRad> & radialData,
		std::vector<double> & coreChargeZ,
		std::vector<VTChg> & frozenCoreElectronChg) const
{
	VTReg buffer;
	Algorithms::helperfunctions::read_binary_file<VTReg, float>(
			filepath_.c_str(), buffer);

	// the first 3 elements are the VASP code version number
	std::int64_t c = 0;
	const int princV = std::floor(buffer[c++]+0.5);
	const int majorV = std::floor(buffer[c++]+0.5);
	const int minorV = std::floor(buffer[c++]+0.5);
	std::cout << "Reading data generated by VASP "<<princV<<"."<<majorV<<"."<<minorV<<"\n";

	// num atoms
	const int numAtoms = std::floor(buffer[c++]+0.5);
	if ( (numAtoms < 1) && (numAtoms > 100000) )
	{
		throw std::runtime_error(std::string("Error processing file ")+
				filepath_.c_str()+": Number of atoms is nonsense (<1 or > 10^5)");
	}

	angularLMaxPerAtom.resize(numAtoms);
	radialPointsPerAtom.resize(numAtoms);
	radiusPerAtom.resize(numAtoms);
	radialData.resize(numAtoms);
	coreChargeZ.resize(numAtoms);
	frozenCoreElectronChg.resize(numAtoms);
	centerAtom.resize(numAtoms);
	namesAtom.resize(numAtoms);
	massAtom.resize(numAtoms);

	// regular grid and its data
	regularGridDim[0] = std::floor(buffer[c++]+0.5);
	regularGridDim[1] = std::floor(buffer[c++]+0.5);
	regularGridDim[2] = std::floor(buffer[c++]+0.5);
	regularData.resize(regularGridDim[0]*regularGridDim[1]*regularGridDim[2]);
	for (int i = 0; i < regularGridDim[0]*regularGridDim[1]*regularGridDim[2]; ++i)
		regularData[i] = buffer[c++];

	for (int atomIndex = 0 ; atomIndex < numAtoms; ++atomIndex)
	{
		// fetch center and atomic mass
		std::array<double,3> center;
		for (auto & xi : center)
			xi = buffer[c++];
		centerAtom[atomIndex] = center;
		massAtom[atomIndex] = buffer[c++];

		// fetch atom symbol
		char firstLatter = static_cast<char>(std::floor(buffer[c++]+0.5));
		char secondLatter = static_cast<char>(std::floor(buffer[c++]+0.5));
		namesAtom[atomIndex] = secondLatter == ' '? std::string({firstLatter}) : std::string({firstLatter, secondLatter});

		// radial grid for this atom
		int numRadPtsGrid = std::floor(buffer[c++]+0.5);
		if(c + 2*numRadPtsGrid + 3 >= buffer.size())
			throw std::runtime_error(std::string("\nError processing file ")+filepath_.c_str()+
							":\n Number of data elements not matching internal description within the file. Probably the file is corrupt.\n");

		radialPointsPerAtom[atomIndex].resize(numRadPtsGrid);
		for (int ir = 0 ; ir < numRadPtsGrid; ++ir)
			radialPointsPerAtom[atomIndex][ir] = buffer[c++];
		radiusPerAtom[atomIndex] = *radialPointsPerAtom[atomIndex].rbegin();

		// angular moment max
		angularLMaxPerAtom[atomIndex] = std::floor(buffer[c++]+0.5);

		// todo implement spin dependence
		const int spinIndex = std::floor(buffer[c++]+0.5);
		if (spinIndex != 1 )
			throw std::runtime_error(std::string("\nError processing file ")+filepath_.c_str()+
					"File shows spin dependence. This functionality is not yet implemented in elephon!");

		// read the core charge, both electronic and core point
		coreChargeZ[atomIndex] = buffer[c++];
		frozenCoreElectronChg[atomIndex].resize(numRadPtsGrid);
		for (int ir = 0 ; ir < numRadPtsGrid; ++ir)
			frozenCoreElectronChg[atomIndex][ir] = buffer[c++];

		// potential data
		const int nAngChnls = std::pow(angularLMaxPerAtom[atomIndex]+1,2);
		if(c + (numRadPtsGrid+2)*nAngChnls > buffer.size())
			throw std::runtime_error(std::string("\nError processing file ")+filepath_.c_str()+
							":\n Number of data elements for the potential are not matching internal description within the file.\n"
							"Probably the file is corrupt.\n");

		radialData[atomIndex].resize(nAngChnls*numRadPtsGrid, 0);
		for (int l = 0 ; l <= angularLMaxPerAtom[atomIndex]; ++l)
			for (int m = -l ; m <= l; ++m)
			{
				const int readL = std::floor(buffer[c++]+0.5);
				const int readM = std::floor(buffer[c++]+0.5);
				const int cnsqChannel = Auxillary::memlayout::angular_momentum_layout(readL, readM);
				for (int ir = 0 ; ir < numRadPtsGrid; ++ir)
				{
					radialData[atomIndex][ir + numRadPtsGrid*cnsqChannel] = buffer[c++];
					if ( std::isnan(std::abs(radialData[atomIndex][ir + numRadPtsGrid*cnsqChannel])) )
						throw std::runtime_error(
								std::string("\nError processing file ")+filepath_.c_str()+":\n" +
								"For atom index "+std::to_string(atomIndex)+" element with angular momentum quantum numbers (l,m)=("+
								std::to_string(l)+","+std::to_string(m)+"), radial index ir="+std::to_string(ir)+" is NaN!\n Exiting ...");
				}
			}
	}
	if (c != buffer.size())
		throw std::runtime_error(std::string("\nError processing file ")+filepath_.c_str()+
			":\n Number of data elements for the potential are not matching internal description within the file.\n"
			"Too many elements loaded for what we can use. Probably the file is corrupt.\n");
}

} /* namespace IOMethods */
} /* namespace elephon */
