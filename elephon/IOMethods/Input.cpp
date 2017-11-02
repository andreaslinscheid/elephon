/*	This file Input.cpp is part of elephon.
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
 *  Created on: Apr 24, 2017
 *      Author: A. Linscheid
 */

#include "IOMethods/Input.h"
#include <stdexcept>
#include <string>
#include <fstream>
#include <iostream>
#include <boost/filesystem.hpp>

namespace elephon
{
namespace IOMethods
{

Input::Input( int argc, char* argv[] )
{
	if ( argc != 2)
		throw std::runtime_error( "\tPass input file as parameter, exiting ... !" );
	inputFile_.read_input_file(argv[1]);

	//Fetch the general configuration variables from the file
	opts_.parse_variables(inputFile_);

	//Flag a warning for all keys that have not been used, as this indicates unintended input
	std::vector<std::string> const unusedKeys = inputFile_.get_list_unread_input_parameters();;
	for ( auto & key : unusedKeys )
		std::cout << "WARNING: Key '"+ key +"' defined in the input file has not been\n"
				" matched to a variable used by the code. Typo?" <<std::endl;

	auto scd = opts_.get_scell();
	if ( scd.size() != 3 )
		throw std::runtime_error(" scell: incorrect size ; must be 3 ");
	if ( (scd[0] < 1) or (scd[1] < 1) or (scd[2] < 1) )
		throw std::runtime_error(" scell: incorrect supercell request ; a dimension cannot be < 1 ");
	//TODO Here could be some consistency checks

	if ( opts_.get_phrange().empty() )
	{
		// determine the phonon range automatically
	}
	else if (opts_.get_phrange().size() == 1)
	{
		if (opts_.get_phrange()[0] <= 0 )
		{
			std::cout << "Problem on input: phrange must be positive number" << std::endl;
			std::exit(0);
		}
	}
	else if (opts_.get_phrange().size() == 2)
	{
		if (opts_.get_phrange()[1] - opts_.get_phrange()[0] <= 0 )
		{
			std::cout << "Problem on input: phrange specify a non-empty, positive range" << std::endl;
			std::exit(0);
		}
	}
	else
	{
		std::cout << "Problem on input: phrange not correctly formated" << std::endl;
		std::exit(0);
	}

	if ( opts_.get_phnpts() < 1 )
	{
		std::cout << "Problem on input: phnpts must be >= 1" << std::endl;
		std::exit(0);
	}

	// convert file paths to absolute ones
	boost::filesystem::path rootdir;
	if ( opts_.get_root_dir().empty() )
		rootdir = boost::filesystem::path(boost::filesystem::current_path());
	else
		rootdir = boost::filesystem::path(opts_.get_root_dir());
	opts_.set_root_dir(rootdir.string());

	auto check_is_relative_path_unix = [] (std::string const & path) {
		if ( not path.empty() )
			return (path.front() != '/');
		return true; // define: empty is a relative path
	};

	boost::filesystem::path elphd_p(opts_.get_elphd());
	if ( check_is_relative_path_unix(opts_.get_elphd()) )
		elphd_p = rootdir / opts_.get_elphd();
	opts_.set_elphd(elphd_p.string());

	// in case this flag is empty, we do not consider an independent dense
	// electronic structure calculation
	if ( not opts_.get_eld().empty() )
	{
		boost::filesystem::path eld_p(opts_.get_eld());
		if ( check_is_relative_path_unix(opts_.get_eld()) )
			eld_p = elphd_p / opts_.get_eld();
		opts_.set_eld(eld_p.string());
	}

	if ( not opts_.get_f_a2F().empty() )
	{
		boost::filesystem::path a2F_p(opts_.get_f_a2F());
		if ( check_is_relative_path_unix(opts_.get_f_a2F()) )
			a2F_p = elphd_p / opts_.get_f_a2F();
		opts_.set_f_a2F(a2F_p.string());
	}

	if ( not opts_.get_f_phdos().empty() )
	{
		boost::filesystem::path phdos_p(opts_.get_f_phdos());
		if ( check_is_relative_path_unix(opts_.get_f_phdos()) )
			phdos_p = elphd_p / opts_.get_f_phdos();
		opts_.set_f_phdos(phdos_p.string());
	}

	if ( not opts_.get_f_kpath().empty() )
	{
		boost::filesystem::path kpath_p(opts_.get_f_kpath());
		if ( check_is_relative_path_unix(opts_.get_f_kpath()) )
			kpath_p = rootdir / opts_.get_f_kpath();
		opts_.set_f_kpath(kpath_p.string());
	}

	if ( not opts_.get_f_bands().empty() )
	{
		if ( not boost::filesystem::exists(opts_.get_f_kpath()))
		{
			std::cout << "Must set f_kpath to valid k path file when"
					" bands calculations are requested!" <<std::endl;
			std::exit(0);
		}
		boost::filesystem::path bands_p(opts_.get_f_bands());
		if ( check_is_relative_path_unix(opts_.get_f_bands()) )
			bands_p = rootdir / opts_.get_f_bands();
		opts_.set_f_bands(bands_p.string());
	}

	if ( not opts_.get_f_ph_bands().empty() )
	{
		if ( not boost::filesystem::exists(opts_.get_f_kpath()))
		{
			std::cout << "Must set f_kpath to valid k path file when"
					" bands calculations are requested!" <<std::endl;
			std::exit(0);
		}
		boost::filesystem::path bands_p(opts_.get_f_ph_bands());
		if ( check_is_relative_path_unix(opts_.get_f_ph_bands()) )
			bands_p = elphd_p / opts_.get_f_ph_bands();
		opts_.set_f_ph_bands(bands_p.string());
	}
}

InputOptions const & Input::get_opts() const
{
	return opts_;
}

} /* namespace IOMethods */
} /* namespace elephon */
