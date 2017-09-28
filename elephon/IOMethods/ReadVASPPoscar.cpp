/*	This file ReadVASPPoscar.cpp is part of elephon.
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
 *  Created on: May 14, 2017
 *      Author: A. Linscheid
 */

#include "ReadVASPPoscar.h"
#include "LatticeStructure/LatticeModule.h"
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <numeric>
#include <assert.h>

namespace elephon
{
namespace IOMethods
{

void ReadVASPPoscar::read_file( std::string filename,
		std::vector<std::pair<std::string, double> > const & atoms  )
{
	std::ifstream file( filename.c_str() );
	if ( ! file.good() )
		throw std::invalid_argument( std::string("File ")+filename+" not readable." );

	//The first line is a comment
	std::string buffer;
	std::getline( file , buffer );

	file.exceptions(std::ifstream::failbit);

	//The second line is 1 or 3 lattice scaling factors
	double scale[3];
	std::getline( file , buffer );
	std::stringstream ssLat(buffer);
	if ( ! (ssLat >> scale[0]) )
		throw std::invalid_argument( std::string("File ")+filename+" has an invalid second line." );
	if ( ssLat >> scale[1]) // true if 3 values are given
	{
		//The other 2 numbers are given too ...
		if ( ! (ssLat >> scale[2]) )
			throw std::invalid_argument( std::string("File ")+filename+" has an invalid second line." );
	}
	else // one value
	{
		//If it is negative VASP treats it as the volume of the cell - we don't accept this feature
		if ( scale[0] < 0 )
			throw std::logic_error( "Cannot handle volume settings in POSCAR, value in second line must be positive" );
		scale[1] = scale[0];
		scale[2] = scale[0];
	}

	//Each of the three lines that follow contains one lattice vector.
	//We store it in our internal lattice matrix A. Note that we use a matrix in the sense
	//that a vector in internal coordinates x becomes a Cartesian vector r by
	// r = A*x
	//So we need to transpose ...
	std::vector<double> latticeMatrix(9);
	for ( int i = 0 ; i < 3 ; ++i )
	{
		std::getline( file , buffer );
		std::stringstream ss(buffer);
		ss.exceptions(std::ifstream::failbit);
		for ( int j = 0 ; j < 3 ; ++j )
		{
			ss >> latticeMatrix[j*3+i];
			latticeMatrix[j*3+i] *= scale[j];
		}
	}

	latticeMatrix_ = std::move(latticeMatrix);

	//Now comes atom type information.
	//This is either a number or a list of atom names
	int nAtomTypes;
	std::vector<std::string> atomsList;
	std::getline( file , buffer );
	std::stringstream ss(buffer);
	if ( ! (ss >> nAtomTypes) )
	{
		std::string next;
		std::stringstream ss(buffer);
		while( ss >> next)
		{
			atomsList.push_back(std::move(next));
		}
	}
	if ( atomsList.size() != 0 )
		nAtomTypes = static_cast<int>(atomsList.size());

	if ( not atomsList.empty() )
	{
		if ( atomsList.size() != atoms.size() )
				throw std::runtime_error("Error parsing POSCAR file: # of atoms and read-in atoms in file different");
		//The list of atoms was provided - we need to read the next line with
		//the number of atoms per atom type which is parsed in the next step
		std::getline( file , buffer );
	}

	//Now that we know how many atom types, read the number of ions per type
	std::vector<int> numAtomPerType(nAtomTypes);
	std::stringstream ssAtoms(buffer);
	ssAtoms.exceptions(std::ifstream::failbit);
	for (int i = 0 ; i < nAtomTypes; ++i)
		ssAtoms >> numAtomPerType[i];

	int totalNumAtoms = std::accumulate(numAtomPerType.begin(),numAtomPerType.end(),0);

	//Here we create a list where for each atom, we tell the type
	std::vector<std::string> typeForAtom(totalNumAtoms);
	int counter = 0;
	for ( int i = 0 ; i < nAtomTypes; ++i)
		for ( int j = 0 ; j < numAtomPerType[i]; ++j)
			typeForAtom[counter++] = atoms[i].first;
	assert(counter==totalNumAtoms);

	//Check for the selective dynamics switch, or the direct/cartesian switch
	bool selectiveDynamics = false;
	bool cartesianInput = false;
	std::getline( file , buffer );
	std::stringstream ssSwitch(buffer);
	std::string word;
	ssSwitch>>word;
	if ( (word.front() == 'S') || (word.front() == 's') )
	{
		selectiveDynamics = true;
		std::getline( file , buffer );
		std::stringstream ssSwitch(buffer);
		ssSwitch>>word;
	}
	if ( (word.front() == 'C') || (word.front() == 'c') || (word.front() == 'k') || (word.front() == 'K') )
		cartesianInput = true;

	//This defines a lambda that can read Fortran bool flags
	//TODO find some reference where a list of allowed Fortran values is given.
	//		I could not find any.
	auto fortran_bool_parser = [] (std::string word)
	{
		if ( (word.compare("T")==0) || (word.compare("t")==0)
				|| (word.compare(".True.")==0) || (word.compare(".true.")==0)
				|| (word.compare(".T.")==0) || (word.compare(".t.")==0) )
			return true;
		return false;
	};

	//Read the coordinates
	for (int i = 0 ; i < totalNumAtoms; ++i)
	{
		std::vector<bool> frozen = {false,false,false};
		std::vector<double> tauI = {0.0,0.0,0.0};
		std::getline( file , buffer );
		std::stringstream ssCoords(buffer);
		ssCoords.exceptions(std::ifstream::failbit);
		for (int xi = 0 ; xi < 3; ++xi)
			ssCoords >> tauI[xi];

		if ( cartesianInput )
		{
			//Apply scaling
			for (int xi = 0 ; xi < 3; ++xi)
				tauI[xi] *= scale[xi];

			LatticeStructure::LatticeModule lattice( latticeMatrix_ );
			lattice.cartesian_to_direct(tauI);
		}

		if ( selectiveDynamics )
		{
			std::string word;
			for (int xi = 0 ; xi < 3; ++xi)
			{
				ssCoords >> word;
				frozen[xi] = fortran_bool_parser(word);
			}
		}
		LatticeStructure::Atom a(atoms[i].second, typeForAtom[i], tauI, frozen);
		atoms_.push_back( std::move(a) );
	}
}

std::vector<LatticeStructure::Atom>
ReadVASPPoscar::get_atoms_list( ) const
{
	return atoms_;
}

std::vector<double> ReadVASPPoscar::get_lattice_matrix() const
{
	return latticeMatrix_;
}

} /* namespace IOMethods */
} /* namespace elephon */
