/*	This file WriteVASPRealSpaceData.cpp is part of elephon.
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
 *  Created on: Jul 2, 2017
 *      Author: A. Linscheid
 */

#include "IOMethods/WriteVASPRealSpaceData.h"
#include <fstream>
#include <ctime>
#include <iomanip>
#include <map>
#include <algorithm>
#include <iostream>
#include <sstream>

namespace elephon
{
namespace IOMethods
{

void
WriteVASPRealSpaceData::write_file(std::string const & filename,
		std::string comment,
		std::vector<int> const & dataDims,
		LatticeStructure::UnitCell const & unitCell,
		std::vector<double> const & data,
		bool spin_resolved,
		bool xmajor) const
{
	assert( data.size() == dataDims[0]*dataDims[1]*dataDims[2]*( spin_resolved ? 2 : 1 ) );

	std::ofstream file( filename.c_str() );
	if ( not file.good() )
		throw std::runtime_error( std::string("Problem opening file ")+ filename);

	if ( comment.empty() )
	{
		auto t = std::time(nullptr);
		comment = std::string("elephon created this file on ") + std::ctime(&t);
	}
	comment.erase(std::remove(comment.begin(), comment.end(), '\n'), comment.end());
	file << comment << '\n';

	//Now comes the unit cell.
	auto floatAccLine = [] (std::vector<double> v, int digits) {
		std::stringstream s;
		s << " "<< std::fixed  << std::setprecision(6) << std::setw(12) << v[0]  << std::setw(12) << v[1]
				<< std::setw(12) << v[2];
		return s.str();
	};

	file << std::fixed << std::setprecision(16) << "   " << unitCell.get_lattice().get_alat() << '\n';
	file << floatAccLine( unitCell.get_lattice().get_lattice_vector(0) , 6 ) << '\n';
	file << floatAccLine( unitCell.get_lattice().get_lattice_vector(1) , 6 ) << '\n';
	file << floatAccLine( unitCell.get_lattice().get_lattice_vector(2) , 6 ) << '\n';

	//List with atom types
	std::string line;
	std::map<std::string,int> ntypes;
	for ( auto a : unitCell.get_atoms_list() )
		ntypes[a.get_kind()]++;
	for ( auto at : ntypes )
		line += at.first + " ";

	file << line + "\n";

	// Number of atoms per type
	line.clear();
	for ( auto a : unitCell.get_atoms_list() )
		if ( ntypes[a.get_kind()] > 0 ) // make sure equal kinds only get listed once
		{
			line += std::to_string(ntypes[a.get_kind()]) + " ";
			ntypes[a.get_kind()] = 0;
		}
	file << line + "\n";

	//Flag for coordinate types
	file << "Direct\n";

	//Positions of atoms
	for ( auto a : unitCell.get_atoms_list() )
	{
		auto pos = a.get_position();
		//map to the cell [0,1[
		for ( auto &pi : pos )
			pi -= std::floor(pi);
		file << std::fixed << std::setprecision(6)  << std::setw(10)
			 << pos[0]  << std::setw(10) << pos[1]  << std::setw(10) << pos[2] << '\n';
	}

	//Empty line
	file << '\n';

	//This lambda shifts the dot one to the left so that we can remove the
	//leading digit. Very unfortunately the C/C++ standard requires the first element of a float to be
	//a digit.
	auto parser_headache = [] (char * num, int length) {
		assert(num[length-1] == '\0');
		std::string exponent;
		exponent.reserve(3);
		int iExp = length-1;
		for ( ; iExp-- ; )
		{
			exponent += num[iExp];
			if ( (num[iExp] == '+') || (num[iExp] == '-') )
				break;
		}
		std::reverse(exponent.begin(), exponent.end() );
		int newExp = std::atoi(exponent.c_str());
		newExp++;

		std::vector<char> newExponent(exponent.size()+2);
		int nwrite = std::snprintf( &newExponent[0], newExponent.size(), "%+03d", newExp);
		newExponent.resize(nwrite);
		if ( newExponent.size() > exponent.size() )
		{
			//We need to delete the last mantissa char to keep the string length.
			char buf[] = {num[iExp-1],num[iExp]};
			assert( (buf[0]=='E') || (buf[0]=='e') );
			num[iExp-2] = buf[0];
			num[iExp-1] = buf[1];
		}
		else if ( newExponent.size() < exponent.size() )
		{
			//We need to append a digit to the mantissa to keep the string length.
			char buf[] = {num[iExp-1],num[iExp]};
			assert( (buf[0]=='E') || (buf[0]=='e') );
			num[iExp-1] = '0';
			num[iExp+0] = buf[0];
			num[iExp+1] = buf[1];
		}

		//Place the new exponent
		for ( int i = 0 ; i < newExponent.size(); ++i )
			num[length-1-newExponent.size()+i] = newExponent[i];

		//Swap the dot and the first number
		for ( int i = 0 ; i < length ; ++i)
			if ( (num[i] >= '0') && (num[i] <= '9') )
			{
				assert( (i+1 < length) && (num[i+1] == '.') );
				std::swap( num[i], num[i+1] );
				break;
			}
	};

	//Header successfully written - write the data
	auto w_data = [&parser_headache] (double const * dataPtr, std::vector<int> const & dataDims) {
		std::string dataLines;
		std::vector<char> buff(512);
		//Fourier dimension
		dataLines += "   "+std::to_string(dataDims[0])+"   "+std::to_string(dataDims[1])+
						"   "+std::to_string(dataDims[2]);
		for (int i = 0 ; i < dataDims[0]*dataDims[1]*dataDims[2]; ++i )
		{
			if ( i % 5 == 0 )
				dataLines += '\n';
			int nWrite = std::snprintf(buff.data(),buff.size(),"%.10E",dataPtr[i]);
			//Annoyingly, c/c++ requires a leading digit such as -0.023 while VASP does not put it (-.023)
			if ( std::isnan(dataPtr[i]) )
				throw std::logic_error("Values to be printed in WriteVASPRealSpaceData::write_file are NaN");
			if ( dataPtr[i] != 0 )
				parser_headache(buff.data(),nWrite+1);
			if ( dataPtr[i] <= 0 )
				dataLines += " "+std::string(buff.data(),buff.data()+nWrite);
			else if ( dataPtr[i] > 0 )
				dataLines += " 0"+std::string(buff.data(),buff.data()+nWrite);
		}
		return dataLines + '\n';
	};

	std::string cntnt;
	if ( xmajor )
	{
		cntnt = w_data( &data[0] , dataDims );
		if ( spin_resolved )
			cntnt += "\n"+w_data( &data[dataDims[0]*dataDims[1]*dataDims[2]] , dataDims );
	}
	else
	{
		std::vector<double> xmajorData(data.size());
		int ns = spin_resolved?2:1;
		for (int is = 0 ; is < ns; ++is)
			for (int iz = 0 ; iz < dataDims[2]; ++iz )
				for (int iy = 0 ; iy < dataDims[1]; ++iy )
					for (int ix = 0 ; ix < dataDims[0]; ++ix )
						xmajorData[((is*dataDims[2]+iz)*dataDims[1]+iy)*dataDims[0]+ix]
								   = data[((is*dataDims[0]+ix)*dataDims[1]+iy)*dataDims[2]+iz];
		cntnt = w_data( &xmajorData[0] , dataDims );
		if ( spin_resolved )
			cntnt += "\n"+w_data( &xmajorData[dataDims[0]*dataDims[1]*dataDims[2]] , dataDims );
	}


	file << cntnt;
}

} /* namespace IOMethods */
} /* namespace elephon */
