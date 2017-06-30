/*	This file ReadVASPxmlFile.cpp is part of elephon.
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
 *  Created on: May 31, 2017
 *      Author: A. Linscheid
 */

#include "ReadVASPxmlFile.h"
#include <boost/foreach.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <algorithm>
#include <iostream>

namespace elephon
{
namespace IOMethods
{

void
ReadVASPxmlFile::parse_file( std::string filename )
{
	if ( filename_.compare(filename) == 0 )
		return;

	filename_ = std::move(filename);
	std::ifstream file( filename_.c_str() );
	if ( ! file.good() )
		throw std::runtime_error( std::string("Problem opening file ")+filename_);

	using boost::property_tree::ptree;
	boost::property_tree::ptree pt;
	read_xml(file, pt);

	std::vector<double> newForces;
	BOOST_FOREACH( ptree::value_type const& val, pt.get_child("modeling.calculation") )
	{
	    if(val.first == "varray")
	    {
	        std::string temp = val.second.get_child("<xmlattr>.name").data();
	        if ( temp == "forces")
	        {
	        	BOOST_FOREACH( ptree::value_type const& val2, val.second )
				{
	        		if ( val2.first == "v" )
	        		{
	        			double tmp;
	        			std::stringstream ss( val2.second.data() );
	        			for ( int i = 0 ; i < 3 ; ++i)
	        			{
	        				ss >> tmp;
	        				newForces.push_back( tmp );
	        			}
	        		}
				}
	        }
	    }

	    if(val.first == "dos")
	    {
	    	BOOST_FOREACH(  ptree::value_type const& val2, val.second )
			{
        		if ( (val2.first == "i") && (val2.second.get_child("<xmlattr>.name").data() == "efermi") )
        		{
        			eFermi_ = strtod( val2.second.data().c_str(), NULL );
        		}

			}
	    }
	};

	forces_ = newForces;

	std::vector<double> newKpts;
	BOOST_FOREACH( ptree::value_type const& val, pt.get_child("modeling.kpoints") )
	{
	    if(val.first == "varray")
	    {
	        std::string temp = val.second.get_child("<xmlattr>.name").data();
	        if ( temp == "kpointlist")
	        {
	        	BOOST_FOREACH( ptree::value_type const& val2, val.second )
				{
					if ( val2.first == "v" )
					{
						double tmp;
						std::stringstream ss( val2.second.data() );
						for ( int i = 0 ; i < 3 ; ++i)
						{
							ss >> tmp;
							newKpts.push_back( tmp );
						}
					}
				}
	        }
	    }
	}
	kpoints_ = newKpts;
}

std::vector<double> const &
ReadVASPxmlFile::get_forces() const
{
	return forces_;
}

std::vector<double> const &
ReadVASPxmlFile::get_k_points() const
{
	return kpoints_;
}

double
ReadVASPxmlFile::get_Fermi_energy() const
{
	return eFermi_;
}

} /* namespace IOMethods */
} /* namespace elephon */
