/*	This file ReadVASPPotcar.cpp is part of elephon.
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

#include <IOMethods/ReadVASPLocpot.h>
#include <fstream>
#include <stdexcept>
#include <sstream>

namespace elephon
{
namespace IOMethods
{

void
ReadVASPLocpot::read_scf_potential(std::string const & filename,
		std::vector<int> & dims,
		std::vector<double> & potential) const
{
	std::ifstream file( filename.c_str() );
	if ( not file.good() )
		throw std::runtime_error( std::string("Problem opening file ")+ filename);

	std::string linebuffer;
	//Comment
	std::getline(file,linebuffer);

	//Now comes the unit cell. In principle this could be compared to input for security
	std::getline(file,linebuffer); // lattice dim
	std::getline(file,linebuffer); //a1
	std::getline(file,linebuffer); //a2
	std::getline(file,linebuffer); //a2

	//List with atom types
	std::getline(file,linebuffer);
	// Number of atoms per type
	std::getline(file,linebuffer);
	std::stringstream ss(linebuffer);
	int numAtoms = 0, nPerType;
	while( ss >> nPerType )
	{
		numAtoms += nPerType;
	}

	if ( numAtoms <= 0 )
		throw std::runtime_error(std::string("Problem parsing file: #atoms not positive in file ")+ filename);

	//Flag for coordinate types
	std::getline(file,linebuffer);

	//skip the positions, too
	for ( int i = 0 ; i < numAtoms ; ++i)
		std::getline(file,linebuffer);
	//Empty line
	std::getline(file,linebuffer);

	//Read fourier dimension
	dims.resize( 3 );
	std::getline(file,linebuffer);
	std::stringstream ssFFT(linebuffer);
	if ( ! (ssFFT >> dims[0] >> dims[1] >> dims[2]) )
		throw std::runtime_error(std::string("Problem parsing file: reading 3 dimensions in file ")+ filename);
	double dummy;
	if ( ssFFT >> dummy)
		throw std::runtime_error(std::string("Problem parsing file: reading more than 3 dimensions in file ")+ filename);

	if ( dims[0]*dims[1]*dims[2] <= 0 )
		throw std::runtime_error(std::string("Problem parsing file: dimensions nonsense in file ")+ filename);

	//Header successfully read - read the data
	std::int64_t beginData = file.tellg();
	file.seekg(0, std::ios::end);
	std::int64_t beginEnd = file.tellg();
	file.seekg(beginData);
	std::int64_t nText = beginEnd-beginData;

	std::vector<char> buffer(nText);
	file.read( buffer.data(), nText );
	std::string content(buffer.begin(), buffer.end());
	std::stringstream dataStream( std::move(content) );

	potential.resize( dims[0]*dims[1]*dims[2] );
	int counter = 0;
	while ( dataStream >> potential[counter] )
	{
		if ( counter == dims[0]*dims[1]*dims[2] )
			throw std::runtime_error(std::string("Problem parsing file: too many data values in file ")+ filename);
		counter++;
	};
	if ( counter != dims[0]*dims[1]*dims[2] )
		throw std::runtime_error(std::string("Problem parsing file: incorrect number of data values in file ")+ filename);
}

} /* namespace IOMethods */
} /* namespace elephon */
