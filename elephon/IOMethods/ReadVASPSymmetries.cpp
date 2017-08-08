/*	This file ReadVASPSymmetries.cpp is part of elephon.
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
 *  Created on: May 16, 2017
 *      Author: A. Linscheid
 */

#include "ReadVASPSymmetries.h"
#include <iterator>
#include <boost/regex.hpp>
#include <fstream>
#include <iostream>
#include <sstream>

namespace elephon
{
namespace IOMethods
{

void ReadVASPSymmetries::read_file(std::string filename )
{
	//Read the OUTCAR content
	std::ifstream file(filename.c_str());
	if ( ! file.good() )
		throw std::runtime_error( std::string("Problem reading file ")+filename);
	file.seekg(0, std::ios::end);
	size_t size = file.tellg();
	std::string filecontent(size, ' ');
	file.seekg(0);
	file.read( &filecontent[0], size);

	//TODO: Find a way to extract that information from the VASP outcar file
    timeReversal_ = true;

	//fetch the blocks with each symmetry
	std::vector<std::string> blocks;
	this->parse_symmetry_blocks( filecontent , blocks );

	if ( blocks.empty() )
	{
		//Check if this file is good, but was run switching symmetries off.
		//VASP does not print the symmetry in this case even though, of cause, there still is
		//idenity operation
		boost::regex nosym(
				"NOSYMM: \\(Re-\\)init\\w* of all symmetry stuff for point group C_1");
		if ( boost::regex_search(filecontent.begin(),filecontent.end(), nosym ) )
		{
		    symmetries_ = std::vector<int>( {1,0,0,0,1,0,0,0,1} );
		    fractionTranslations_ = std::vector<double>( {0.0, 0.0, 0.0} );
		    timeReversal_ = false;
		    return;
		}
		else
			throw std::runtime_error( std::string("Unable to parse symmetries from file ")+filename );
	}

	const char * reSym =
			"\\s*isymop:\\s*([+-]?\\d+)\\s+([+-]?\\d+)\\s+([+-]?\\d+)\\s*\\n"
			"\\s*\\s*([+-]?\\d+)\\s+([+-]?\\d+)\\s+([+-]?\\d+)\\s*\\n"
			"\\s*\\s*([+-]?\\d+)\\s+([+-]?\\d+)\\s+([+-]?\\d+)\\s*\\n";
    boost::regex symmetryRegex(reSym);

	const char * reFract =
			"\\s*gtrans:\\s*([+-]?\\d+\\.\\d+)\\s+([+-]?\\d+\\.\\d+)\\s+([+-]?\\d+\\.\\d+)\\s*\\n";
    boost::regex fractionalTransRegex(reFract);

    symmetries_ = std::vector<int>( blocks.size()*9 );
    fractionTranslations_ = std::vector<double>( blocks.size()*3 );
	for ( int isym = 0 ; isym < static_cast<int>(blocks.size()); ++isym)
	{
		std::vector< std::string > symmetriesStr;
		boost::match_results<std::string::const_iterator> res;
		boost::regex_search(blocks[isym], res, symmetryRegex );
	    assert( res.size() == 9+1 );
	    for ( int i = 0 ; i < 3 ; ++i)
		    for ( int j = 0 ; j < 3 ; ++j)
		    	symmetries_[(isym*3+i)*3+j] = std::stoi( std::string(res[i*3+j+1].first,res[i*3+j+1].second) );

		std::vector< std::string > fractStr;
		boost::regex_search(blocks[isym], res, fractionalTransRegex );
	    assert( res.size() == 3+1 );
	    for ( int i = 1 ; i < 4 ; ++i)
	    	fractionTranslations_[isym*3+i-1] = std::stod(  std::string(res[i].first,res[i].second)  ) ;
	}
}

void ReadVASPSymmetries::parse_symmetry_blocks(std::string const & fcontent,
		std::vector<std::string> & blocks) const
{
	//TODO this regular expression has likely to be hardened against different versions of VASP
	const char * re =
			//We search the file for the following sequence indicating the symmetries block
			"(irot\\s*\\:\\s*\\d+\\s*\\n"
			"\\s*--------------------------------------------------------------------\\s*\\n"
			"\\s*isymop:\\s*[+-]?\\d+\\s+[+-]?\\d+\\s+[+-]?\\d+\\s*\\n"
			"\\s*\\s*[+-]?\\d+\\s+[+-]?\\d+\\s+[+-]?\\d+\\s*\\n"
			"\\s*\\s*[+-]?\\d+\\s+[+-]?\\d+\\s+[+-]?\\d+\\s*\\n"
			"\\s*\\n"
			"\\s*gtrans:\\s*[+-]?\\d+\\.\\d+\\s+[+-]?\\d+\\.\\d+\\s+[+-]?\\d+\\.\\d+\\s*\\n)";
    boost::regex expression(re);

    std::copy(boost::sregex_token_iterator(fcontent.begin(), fcontent.end(), expression),
        boost::sregex_token_iterator(),
        std::back_inserter(blocks));
}

std::vector<int> const&
ReadVASPSymmetries::get_symmetries() const
{
	return symmetries_;
}

std::vector<double> const&
ReadVASPSymmetries::get_fractionTranslations() const
{
	return fractionTranslations_;
}

bool
ReadVASPSymmetries::get_time_revesal_symmetry() const
{
	return timeReversal_;
}

} /* namespace IOMethods */
} /* namespace elephon */
