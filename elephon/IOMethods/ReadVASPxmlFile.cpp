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
#include <boost/algorithm/string.hpp>
#include <boost/property_tree/xml_parser.hpp>
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
	read_xml(file, pt_);
	parseForces_ = true;
	parseAtoms_ = true;
	parseEnergies_ = true;
	parseKPoints_ = true;
	parseFFTDims_ = true;
}

std::vector<double> const &
ReadVASPxmlFile::get_forces()
{
	if ( parseForces_ )
		this->parse_forces();
	return forces_;
}

void
ReadVASPxmlFile::parse_forces()
{
	using boost::property_tree::ptree;
	auto calculationLeaf = pt_.get_child("modeling.calculation");
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
	};

	if (newForces.empty())
		throw std::runtime_error(std::string("Error parsing")+filename_+": Unable to find the atomic force data");
	forces_ = newForces;
	parseForces_ = false;
}

void
ReadVASPxmlFile::parse_latticeMat()
{
	using boost::property_tree::ptree;
	std::vector<double> lmat;
	BOOST_FOREACH( ptree::value_type const& val0, pt_.get_child("modeling") )
	{
		if(val0.first == "structure")
		{
	        std::string temp = val0.second.get_child("<xmlattr>.name").data();
	        if ( temp == "finalpos")
	        {
				BOOST_FOREACH( ptree::value_type const& val, val0.second.get_child("crystal") )
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
	        }
		}
	}
	latticeMat_.resize(9);
	for ( int  i = 0 ; i < 3 ; ++i)
		for ( int  j = 0 ; j < 3 ; ++j)
			latticeMat_[i*3+j] = lmat[j*3+i];
	parseLatticeMatrix_ = false;
}

void
ReadVASPxmlFile::parse_atoms()
{
	using boost::property_tree::ptree;

	std::vector<std::vector<double>> atomicPos;
	BOOST_FOREACH( ptree::value_type const& val, pt_.get_child("modeling.calculation.structure") )
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
	BOOST_FOREACH( ptree::value_type const& val, pt_.get_child("modeling.atominfo") )
	{
	    if(val.first == "array")
	    	if(val.second.get_child("<xmlattr>.name").data().compare("atomtypes") == 0 )
	    	{
	    		int c = 0;
				BOOST_FOREACH( ptree::value_type const& val2, val.second.get_child("set") )
				{
					if ( val2.first == "rc" )
					{
						c++; // vasp counts from 1
						BOOST_FOREACH( ptree::value_type const& val3, val2.second )
						{
							if ( val3.first == "c" )
							{
								std::string atomName = val3.second.data();
								boost::algorithm::trim(atomName);
								typeInfo[c].push_back(atomName);
							}
						}
						if ( typeInfo[c].size() < 3 )
							throw std::runtime_error("Problem parsing atom info");
					}
				}
	    	}
	}

	atoms_.clear();
	BOOST_FOREACH( ptree::value_type const& val, pt_.get_child("modeling.atominfo") )
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
							{
								std::string atomName = val3.second.data();
								boost::algorithm::trim(atomName);
								atomInfo.push_back(atomName);
							}
						}
						if (atomInfo.size() != 2)
							throw std::runtime_error("Problem parsing atominfo from vasprun.xml");
						int numAtom = atoms_.size();
						int typenum = std::stoi(atomInfo[1]);
						auto it = typeInfo.find(typenum);
						if ( it == typeInfo.end())
							throw std::runtime_error("Problem finding atom typeinfo.");
						double mass = std::stof(it->second.at(2));
						if ((numAtom<0)||(numAtom>= atomicPos.size()))
							throw std::runtime_error("Problem with atom number outside of range of atomic positions.");
						atoms_.push_back( LatticeStructure::Atom(mass, atomInfo[0], atomicPos[numAtom], {false, false, false}, 1e-6) );
					}
				}
	    	}
	}
	parseAtoms_ = false;
}

std::vector<double> const &
ReadVASPxmlFile::get_k_points()
{
	if (parseKPoints_)
		this->parse_kpoints();
	return kpoints_;
}

double
ReadVASPxmlFile::get_Fermi_energy()
{
	if (parseEnergies_)
		this->parse_energies();
	return eFermi_;
}

std::vector<double> const &
ReadVASPxmlFile::get_energies()
{
	if (parseEnergies_)
		this->parse_energies();
	return energies_;
}

int
ReadVASPxmlFile::get_nBnd()
{
	if (parseEnergies_)
		this->parse_energies();
	return nBnd_;
}

int
ReadVASPxmlFile::get_nkp()
{
	if (parseKPoints_)
		this->parse_kpoints();
	return kpoints_.size()/3;
}

std::vector<double> const &
ReadVASPxmlFile::get_lattice_matrix()
{
	if (parseLatticeMatrix_)
		this->parse_latticeMat();
	return latticeMat_;
}

std::vector<LatticeStructure::Atom> const &
ReadVASPxmlFile::get_atoms_list()
{
	if (parseAtoms_)
		this->parse_atoms();
	return atoms_;
}

std::vector<int>
ReadVASPxmlFile::get_wfct_fourier_dim()
{
	if (parseFFTDims_)
		this->parse_fftdims();
	return wfctFourierDim_;
}

std::vector<int>
ReadVASPxmlFile::get_charge_fourier_dim()
{
	if (parseFFTDims_)
		this->parse_fftdims();
	return chargeFourierDim_;
}

std::vector<int>
ReadVASPxmlFile::get_k_grid_dim()
{
	if (parseKPoints_)
		this->parse_kpoints();
	return kDim_;
}

std::vector<double>
ReadVASPxmlFile::get_k_grid_shift()
{
	if (parseKPoints_)
		this->parse_kpoints();
	return kShift_;
}

void
ReadVASPxmlFile::parse_energies()
{
	using boost::property_tree::ptree;
	auto calculationLeaf = pt_.get_child("modeling");
	BOOST_FOREACH( ptree::value_type const& val, pt_.get_child("modeling") )
	{
	    if(val.first == "calculation")
	    {
	    	boost::optional< const ptree& > child = val.second.get_child_optional("eigenvalues");
	    	if ( child )
	    		calculationLeaf = val.second;
	    }
	}

	BOOST_FOREACH( ptree::value_type const& val, calculationLeaf )
	{
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

	std::vector<double> bandBuffer;
	energies_.clear();
	nBnd_ = 0;

	// We need the k points to make sense of the energies
	this->parse_kpoints();

	int nkp = kpoints_.size()/3;
	BOOST_FOREACH( ptree::value_type const& val, calculationLeaf.get_child("eigenvalues.array.set") )
	{
	    if(val.first == "set")
	    {
	        std::string temp = val.second.get_child("<xmlattr>.comment").data();
			if ( temp == "spin 1")
			{
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
	parseEnergies_ = false;
}

void
ReadVASPxmlFile::parse_kpoints()
{
	using boost::property_tree::ptree;
	std::vector<double> newKpts;
	BOOST_FOREACH( ptree::value_type const& val, pt_.get_child("modeling.kpoints") )
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

	kDim_.assign(3, -1);
	kShift_.assign(3, 0.0);
	BOOST_FOREACH( ptree::value_type const& val, pt_.get_child("modeling.kpoints.generation") )
	{
		if (  val.first != "v")
			continue;
        std::string temp = val.second.get_child("<xmlattr>.name").data();
        if ( temp == "divisions" )
        {
        	std::stringstream ss(val.second.data());
        	ss >> kDim_[0] >> kDim_[1] >> kDim_[2];
        }
        if ( (temp == "usershift") or (temp == "shift") )
        {
        	std::stringstream ss(val.second.data());
        	std::vector<double> kShiftLocal(3,0);
        	ss >> kShiftLocal[0] >> kShiftLocal[1] >> kShiftLocal[2];
        	for ( int i = 0 ; i < 3 ; ++i )
        		kShift_[i] += kShiftLocal[i];
        }
	}
	if ( *std::min_element(kDim_.begin(), kDim_.end()) <= 0 )
		throw std::runtime_error("Problem parsing k grid dimensions from vasprun.xml");

	parseKPoints_ = false;
}

void
ReadVASPxmlFile::parse_fftdims()
{
	using boost::property_tree::ptree;
	wfctFourierDim_.assign(3, 0);
	chargeFourierDim_.assign(3, 0);
	BOOST_FOREACH( ptree::value_type const& val, pt_.get_child("modeling.parameters") )
	{
	    if(val.first == "separator")
	    {
	        std::string temp = val.second.get_child("<xmlattr>.name").data();
			if ( temp == "grids")
			{
				BOOST_FOREACH( ptree::value_type const& val2, val.second )
				{
					if ( val2.first == "i" )
					{
						std::string gridtag = val2.second.get_child("<xmlattr>.name").data();
						if ( gridtag == "NGX" )
							wfctFourierDim_[0] = std::stoi(val2.second.data());
						if ( gridtag == "NGY" )
							wfctFourierDim_[1] = std::stoi(val2.second.data());
						if ( gridtag == "NGZ" )
							wfctFourierDim_[2] = std::stoi(val2.second.data());
						if ( gridtag == "NGXF" )
							chargeFourierDim_[0] = std::stoi(val2.second.data());
						if ( gridtag == "NGYF" )
							chargeFourierDim_[1] = std::stoi(val2.second.data());
						if ( gridtag == "NGZF" )
							chargeFourierDim_[2] = std::stoi(val2.second.data());
					}
				}
			}
	    }
	}
	if ( *std::min_element(wfctFourierDim_.begin(), wfctFourierDim_.end()) <= 0 )
		throw std::runtime_error("Problem parsing fourier grids for wavefunctions from vasprun.xml");
	if ( *std::min_element(chargeFourierDim_.begin(), chargeFourierDim_.end()) <= 0 )
		throw std::runtime_error("Problem parsing fourier grids for the charge from vasprun.xml");

}

} /* namespace IOMethods */
} /* namespace elephon */
