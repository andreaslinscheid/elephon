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

#include <IOMethods/Input.h>
#include <stdexcept>
#include <string>
#include <fstream>
#include <boost/program_options.hpp>

namespace elephon {
namespace IOMethods {

Input::Input(int argc, char const * const* argv)
{
	namespace po = boost::program_options;
	po::options_description desc("Options");
	std::string buf;
	desc.add_options()
		("scellx", po::value<size_t >()->default_value(1),
					"Number of cells in direction x. Default: 1" )
		("scelly", po::value<size_t >()->default_value(1),
					"Number of cells in direction y. Default: 1" )
		("scellz", po::value<size_t >()->default_value(1),
					"Number of cells in direction z. Default: 1" )
		("numFS", po::value<size_t >()->default_value(1000),
					"Number points sampling the Fermi surface. Default: 1000" )
		("wavef", po::value<std::string>()->default_value("WAVECAR"),
					"File with wavefunction data. WAVECAR." )
	;

	if ( argc != 1 )
		throw std::invalid_argument( "Need one parameter as input file" );

	std::ifstream ifs( argv[1] );
	if (!ifs)
		throw std::invalid_argument( std::string("Cannot open input file: ")+ argv[1] );
	po::store(po::parse_config_file(ifs, desc), vm_);
	notify(vm_);
}


std::vector<size_t> Input::get_super_cell_size() const
{
	if (vm_.count("scellx")+vm_.count("scelly")+vm_.count("scellz")<3)
		throw std::logic_error("Supercell parameters must be present");
	return std::vector<size_t>({vm_["scellx"].as<size_t>(),
								vm_["scelly"].as<size_t>(),
								vm_["scellz"].as<size_t>()});
}

size_t Input::get_numFS() const
{
	if (not vm_.count("numFS"))
		throw std::logic_error("Supercell parameters must be present");
	return vm_["numFS"].as<size_t>();
}

} /* namespace IOMethods */
} /* namespace elephon */
