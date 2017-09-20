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
#include <boost/optional/optional.hpp>
#include <algorithm>
#include <iostream>
#include <map>

namespace elephon
{
namespace IOMethods
{

bool
ReadVASPxmlFile::is_parsed() const
{
	return (! filename_.empty());
}

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

	// at this point we determine the last modeling.calculation node, where the data
	// will be read from
	auto calculationLeaf = pt.get_child("modeling");
	BOOST_FOREACH( ptree::value_type const& val, pt.get_child("modeling") )
	{
	    if(val.first == "calculation")
	    {
	    	boost::optional< const ptree& > child = val.second.get_child_optional("eigenvalues");
	    	if ( child )
	    		calculationLeaf = val.second;
	    }
	}
	if ( calculationLeaf == pt.get_child("modeling"))
		throw std::runtime_error("Could not locate the caluclation with the eigenvalues in vasprun.xml");

	std::vector<double> lmat;
	BOOST_FOREACH( ptree::value_type const& val, calculationLeaf.get_child("structure.crystal") )
	{
	    if(val.first == "varray")
	    {
	        std::string temp = val.second.get_child("<xmlattr>.name").data();
	        if ( temp == "basis")
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
	        				lmat.push_back( tmp );
	        			}
	        		}
				}
	        	if ( lmat.size() != 9 )
	        		throw std::runtime_error("Problem reading lattice matrix - size incorrect");
	        }
	    }
	}
	latticeMat_.resize(9);
	for ( int  i = 0 ; i < 3 ; ++i)
		for ( int  j = 0 ; j < 3 ; ++j)
			latticeMat_[i*3+j] = lmat[j*3+i];

	std::vector<std::vector<double>> atomicPos;
	BOOST_FOREACH( ptree::value_type const& val, calculationLeaf.get_child("structure") )
	{
	    if(val.first == "varray")
	    {
	        std::string temp = val.second.get_child("<xmlattr>.name").data();
			if ( temp == "positions")
			{
				BOOST_FOREACH( ptree::value_type const& val2, val.second )
				{
					if ( val2.first == "v" )
					{
						std::vector<double> p(3);
						std::stringstream ss( val2.second.data() );
						for ( int i = 0 ; i < 3 ; ++i)
							ss >> p[i];
						atomicPos.push_back( std::move(p) );
					}
				}
			}
	    }
	}
	std::map<int,std::vector<std::string>> typeInfo;
	BOOST_FOREACH( ptree::value_type const& val, pt.get_child("modeling.atominfo") )
	{
	    if(val.first == "array")
	    	if(val.second.get_child("<xmlattr>.name").data().compare("atomtypes") == 0 )
	    	{
	    		int c = 0;
				BOOST_FOREACH( ptree::value_type const& val2, val.second.get_child("set") )
				{
					if ( val2.first == "rc" )
					{
						c++;
						BOOST_FOREACH( ptree::value_type const& val3, val2.second )
						{
							if ( val3.first == "c" )
								typeInfo[c].push_back(val3.second.data());
						}
					}
				}
	    	}
	}
	atoms_.clear();
	BOOST_FOREACH( ptree::value_type const& val, pt.get_child("modeling.atominfo") )
	{
	    if(val.first == "array")
	    	if(val.second.get_child("<xmlattr>.name").data().compare("atoms") == 0 )
	    	{
				BOOST_FOREACH( ptree::value_type const& val2, val.second.get_child("set") )
				{
					if ( val2.first == "rc" )
					{
			    		std::vector<std::string> atomInfo;
						BOOST_FOREACH( ptree::value_type const& val3, val2.second )
						{
							if ( val3.first == "c" )
								atomInfo.push_back(val3.second.data());
						}
						if (atomInfo.size() != 2)
							throw std::runtime_error("Problem parsing atominfo from vasprun.xml");
						int type = std::stoi(atomInfo[1]);
						int numAtom = atoms_.size();
						atoms_.push_back( LatticeStructure::Atom(atomInfo[0], atomicPos[numAtom], {false, false, false}, 1e-6) );
					}
				}
	    	}
	}

	std::vector<double> newForces;
	std::vector<double> newEV;
	BOOST_FOREACH( ptree::value_type const& val, calculationLeaf )
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
							//apply the elephon 1.BZ convention [-0.5,0.5[. NOTE: VASP uses ]-0.5,0.5]
							tmp -= std::floor(tmp+0.5);
							newKpts.push_back( tmp );
						}
					}
				}
	        }
	    }
	}
	kpoints_ = newKpts;

	std::vector<double> bandBuffer;
	energies_.clear();
	nBnd_ = 0;
	int nkp = kpoints_.size()/3;
	BOOST_FOREACH( ptree::value_type const& val, calculationLeaf.get_child("eigenvalues.array.set") )
	{
	    if(val.first == "set")
	    {
	        std::string temp = val.second.get_child("<xmlattr>.comment").data();
			if ( temp == "spin 1")
			{
				int ispin = 0;
				BOOST_FOREACH( ptree::value_type const& val2, val.second )
				{
					if ( val2.first == "set" )
					{
						std::string ktag = val2.second.get_child("<xmlattr>.comment").data();
						ktag.erase(0,6); // remove the first label 'kpoint'
						int ik = std::stoi(ktag)-1;
						if ( (ik < 0) or (ik >= nkp) )
							throw std::runtime_error("Problem parsing eigenvalues from vasprun.xml - kp out of range");
						int ib = 0;
						BOOST_FOREACH( ptree::value_type const& val3, val2.second )
						{
							if ( val3.first == "r" )
							{
								if ( nBnd_ == 0 )
									bandBuffer.push_back( std::stod(val3.second.data()) );
								else
									energies_[ik*nBnd_ + ib] = std::stod(val3.second.data());
								++ib;
							}
						}
						if ( energies_.empty() )
						{
							nBnd_ = ib;
							energies_.assign(nBnd_*nkp, std::numeric_limits<double>::infinity());
							for ( ib = 0 ; ib < nBnd_ ; ++ib)
								energies_[ik*nBnd_ + ib] = bandBuffer[ib];
						}
					}
				}
			}
	    }
	}
	auto m = std::max_element(energies_.begin(), energies_.end());
	if ( *m == std::numeric_limits<double>::infinity() )
		throw std::runtime_error("Problem parsing energy values from vasprun.xml - energies missing");
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

std::vector<double> const &
ReadVASPxmlFile::get_energies() const
{
	return energies_;
}

int
ReadVASPxmlFile::get_nBnd() const
{
	return nBnd_;
}

int
ReadVASPxmlFile::get_nkp() const
{
	return kpoints_.size()/3;
}

std::vector<double> const &
ReadVASPxmlFile::get_lattice_matrix() const
{
	return latticeMat_;
}

std::vector<LatticeStructure::Atom> const &
ReadVASPxmlFile::get_atoms_list() const
{
	return atoms_;
}

} /* namespace IOMethods */
} /* namespace elephon */
