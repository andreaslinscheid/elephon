/*	This file KPath.cpp is part of elephon.
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
 *  Created on: Nov 1, 2017
 *      Author: A. Linscheid
 */

#include "IOMethods/KPath.h"
#include <fstream>
#include <stdexcept>
#include <cmath>
#include <sstream>
#include <boost/regex.hpp>
#include <algorithm>
#include <assert.h>

namespace elephon
{
namespace IOMethods
{

void
KPath::read_kpath_file(std::string const & filename)
{
	std::ifstream file( filename.c_str() );
	if ( not file.good() )
		throw std::runtime_error(std::string()+"Problem opening kpath file "+filename);

	//figure out the dimension of the input k points
	std::string firstLine;
	std::getline(file,firstLine);
	std::stringstream ss(firstLine);
	std::string dummy;
	int nelem = 0;
	while( ss >> dummy)
	{
		nelem++;
	}
	int dim = (nelem-2-1)/2;
	if ( dim != 3 )
		throw std::runtime_error("Problem reading file: dimension of k space figured out !=3");
	file.clear();
	file.seekg(0,std::ios::beg);

	std::string label1, label2;
	std::vector<double> v1(3,0);
	std::vector<double> v2(3,0);
	std::vector<double> kpath;
	while( true )
	{
		//read first k point
		if (not ( file >> label1 ) )
		{
			if ( file.eof() )
				break;
			throw std::runtime_error(std::string()+"Problem reading label1 kpath file "+filename);
		}

		for ( int i = 0 ; i < 3; i++)
			if (not ( file >> v1[i] ) )
				throw std::runtime_error(std::string()+"Problem reading vector1 kpath file "+filename);

		int nkpts_sect = 0;
		if (not ( file >> nkpts_sect ) )
			throw std::runtime_error(std::string()+"Problem reading section nkpts kpath file "+filename);

		//read second k point
		if (not ( file >> label2 ) )
			throw std::runtime_error(std::string()+"Problem reading label2 kpath file "+filename);

		for ( size_t i = 0 ; i < 3; i++)
			if (not ( file >> v2[i] ) )
				throw std::runtime_error(std::string()+"Problem reading vector2 kpath file "+filename);

		labels_.push_back( std::make_pair(kpath.size()/3,label1) );
		for ( int iks = 0 ; iks < nkpts_sect; iks++)
			for ( int i = 0 ; i < 3; i++)
			{
				double vfbz = v1[i]+(v2[i]-v1[i])*static_cast<double>(iks)/static_cast<double>(nkpts_sect);
				vfbz -= std::floor(vfbz + 0.5);
				kpath.push_back( vfbz );
			}
	}

	//add the last k points
	labels_.push_back( std::make_pair(kpath.size()/3,label2) );
	for ( int i = 0 ; i < 3; i++)
	{
		double vfbz = v2[i];
		vfbz -= std::floor(vfbz + 0.5);
		kpath.push_back( vfbz );
	}
	kpts_.assign(kpath.begin(),kpath.end());
}

void
KPath::produce_gnuplot_script_spectral_function(
		std::string filenameScript,
		std::string filenameData,
		std::string quantityLabel,
		std::vector<double> const & frequencies,
		std::vector<double> const & data) const
{
	auto mm = std::minmax_element(frequencies.begin(), frequencies.end());
	std::string plotscript =
			this->build_plot_script_common(quantityLabel,
					std::make_pair(*mm.first, *mm.second)) +
	"unset colorbox\n"
	"set pm3d map\n"
	"set palette model RGB functions r(gray), g(gray), b(gray)\n"
	"cutOffset(x) = x > 0.0 ? x : 0.0\n"
	"r(x)=cutOffset(-2.5*x+0)   + 0.8*cutOffset(-1+2.5*abs(x))\n"
	"g(x)=cutOffset(-2.5*x+1)   + 0.8*cutOffset(-1+2.5*abs(x))\n"
	"b(x)=cutOffset(-2.5*x+2)   + 0.8*cutOffset(-1+2.5*abs(x))\n"
	"sp '"+filenameData+"' binary matrix using 2:1:(abs($3)>exp(0)?0:-log(abs($3))>10?10:-log(abs($3))) w pm3d";

	this->dump_file(filenameScript, plotscript);
}

void
KPath::produce_gnuplot_script_stable_particle(
		std::string filenameScript,
		std::string filenameData,
		std::string quantityLabel,
		std::vector<double> const & data,
		int numBnds,
		std::pair<double,double> energyRange) const
{
	this->save_bands_path_stable_particle(filenameData, data, numBnds);
	std::string plotscript = this->build_plot_script_common(quantityLabel, energyRange);
	plotscript += "unset key\n";
	std::string plotline =
			"p '"+filenameData+"'";
	for (int ib = 0 ; ib < numBnds ; ++ib)
		plotline += " u 1:"+std::to_string(ib+2)+" w l lt 1 lc rgb \"black\", ''";
	plotline.erase(static_cast<int>(plotline.length())-4);
	plotscript += plotline + "\n";

	this->dump_file(filenameScript, plotscript);
}

void
KPath::save_bands_path_stable_particle(
		std::string const & filename,
		std::vector<double> const & data,
		int numBnds) const
{
	int numK = data.size()/numBnds;
	assert(numK == kpts_.size()/3);
	assert(kpts_.size()/3>1);
	std::string fileContent = "#bands units are eV";
	for (int ikpath = 0 ; ikpath < numK ; ++ikpath)
	{
		fileContent += std::to_string(static_cast<double>(ikpath)/static_cast<double>(numK));
		for (int ib = 0 ; ib < numBnds ; ++ib)
			fileContent += "\t"+std::to_string(data[ikpath*numBnds+ib]);
		fileContent += "\n";
	}
	this->dump_file(filename, fileContent);
}

std::vector<double> const &
KPath::get_k_points() const
{
	return kpts_;
}

std::string
KPath::build_plot_script_common(
		std::string quantityLabel,
		std::pair<double,double> minMax) const
{
	std::string xTics = this->build_xtics();

	quantityLabel = boost::regex_replace(quantityLabel, boost::regex("\\\\"), "\\\\\\\\");

	std::string plotscript =
	"#!/usr/bin/gnuplot\n"
	"reset\n"
	"set terminal epslatex color font 'Helvetica,10' lw 2\n"
	"set output './tmp.tex'\n\n"
	"unset key\n"
	"set grid xtics\n"
	"set xtics ("+xTics+")\n\n"
	"set ylabel \""+quantityLabel+"\" offset 0,0,0\n"
	"unset xlabel\n\n"
	"set yrange ["+std::to_string(minMax.first)+":"+std::to_string(minMax.second)+"]\n\n"
	"set arrow from 0,0 to 1,0 lw 1 lt 2 lc rgb \"black\" nohead front\n"; // Fermi level
	return plotscript;
}

std::string
KPath::build_xtics() const
{
	//build the path label
	std::string xTics;
	boost::regex escape("\\\\");
	for ( auto l : labels_)
	{
		double pos = (double(l.first))/double(kpts_.size()/3);
		// gnuplot skips values on the very end of the range 0 - 1 ... annoying.
		pos = std::max(pos, 0.01);
		pos = std::min(pos, 0.99);
		l.second = boost::regex_replace(l.second, escape, "\\\\\\\\"); // escape sequence ...
		xTics += std::string("\"") + l.second
				+"\" "+std::to_string(pos)+",";
	}
	xTics.pop_back();
	return xTics;
}

void
KPath::dump_file(
		std::string const & filename,
		std::string const & content) const
{
	std::ofstream scriptFile(filename.c_str());
	if ( ! scriptFile.good() )
		throw std::runtime_error( std::string("Error creating file ")+filename);
	scriptFile << content;
	scriptFile.close();
}

} /* namespace IOMethods */
} /* namespace elephon */
