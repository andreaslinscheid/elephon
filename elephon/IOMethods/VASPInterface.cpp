/*	This file VASPInterface.cpp is part of elephon.
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
 *  Created on: May 17, 2017
 *      Author: A. Linscheid
 */

#include "VASPInterface.h"
#include "IOMethods/InputFile.h"
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <map>
#include <sstream>
#include <iomanip>
#include <fstream>

namespace elephon
{
namespace IOMethods
{

void VASPInterface::set_up_run(
		std::string root_directory,
		std::string target_directory,
		std::vector<int> const & kptSampling,
		std::vector<double> const & kptShift,
		LatticeStructure::UnitCell const & unitcell,
		std::map<std::string,std::string> const& options) const
{
	boost::filesystem::path root(root_directory);
	boost::filesystem::path elphd(target_directory);

	std::vector<std::string> baseFiles = this->list_all_input_files();

	//first things first, we copy the POTCAR file to the location
	boost::filesystem::path potcarPrev = root / baseFiles[3] ;
	boost::filesystem::path potcarNew = elphd /  baseFiles[3];
	boost::filesystem::copy( potcarPrev, potcarNew );

	//write parameters in POSCAR according to data in unitcell
	auto atomOrder = this->read_potcar_atom_order( potcarNew.string() );
	boost::filesystem::path poscarNew = elphd / baseFiles[1];
	this->overwrite_POSCAR_file( poscarNew.string(), atomOrder, unitcell );

	//Write the KPOINTS file
	boost::filesystem::path kpts = elphd / baseFiles[2];
	this->write_KPOINTS_file( kpts.string(), kptShift, kptSampling );

	//First copy and then modify the INCAR file
	boost::filesystem::path incarPrev = root / baseFiles[0] ;
	boost::filesystem::path incarNew = elphd /  baseFiles[0];
	boost::filesystem::copy( incarPrev, incarNew );

	this->modify_incar_file( incarNew.string(), options );
}

std::map<std::string,std::string>
VASPInterface::options_nscf_keep_wfctns_no_relax() const
{
	std::map<std::string,std::string> options;
	//We want to keep the wavefunctions
	options["LWAVE"] = ".TRUE.";
	options["IBRON"] = "-1";
	options["NSW"] = "0";
	options["ICHARG"] = "11";
	return options;
}

std::map<std::string,std::string>
VASPInterface::options_scf_supercell_no_wfctns_no_relax() const
{
	std::map<std::string,std::string> options;
	//We want to keep the wavefunctions
	options["LWAVE"] = ".FALSE.";
	options["IBRON"] = "-1";
	options["NSW"] = "0";
	options["ICHARG"] = "1";
	options["LVTOT"] = ".TRUE.";
	options["PREC"] = "High";
	options["LVHAR"] = ".FALSE.";
	return options;
}

void VASPInterface::read_wavefunctions(
		std::vector<std::string> const & files,
		std::vector<int> const & kpts,
		std::vector<int> const & bandIndices,
		std::vector< std::complex<float> > & wfctData,
		std::vector< std::vector<int> > & fourierMap,
		std::vector<int> & fftDim)
{
	wfcReader_.prepare_wavecar( files.at(0) );
	wfcReader_.read_wavefunction( kpts, bandIndices, wfctData, fourierMap, fftDim );
}

void VASPInterface::read_atoms_list(
		std::vector<std::string> const & baseFiles,
		std::vector<LatticeStructure::Atom> & atoms)
{
	if ( posReader_.get_atoms_list().empty() )
		posReader_.read_file( baseFiles.at(0) );

	atoms = posReader_.get_atoms_list();
}

void VASPInterface::read_cell_paramters(
		std::vector<std::string> const & baseFiles,
		LatticeStructure::LatticeModule & lattice)
{
	if ( posReader_.get_atoms_list().empty() )
		posReader_.read_file( baseFiles.at(0) );

	lattice.initialize( posReader_.get_lattice_matrix() );
}

void VASPInterface::read_electronic_potential(
		std::vector<std::string> const & files,
		std::vector<float> & output)
{

}

void VASPInterface::read_symmetries(
		std::vector<std::string> const & files,
		double symPrec,
		LatticeStructure::Symmetry & symmetry)
{
	if ( symReader_.get_symmetries().empty() )
		symReader_.read_file( files.at(0) );

	symmetry.initialize(symPrec, symReader_.get_symmetries(), symReader_.get_fractionTranslations() );
}

std::vector<std::string> VASPInterface::list_all_input_files() const
{
	std::vector<std::string> list;
	list.push_back("INCAR");
	list.push_back("POSCAR");
	list.push_back("KPOINTS");
	list.push_back("POTCAR");
	return list;
}

std::vector<std::string >
VASPInterface::read_potcar_atom_order( std::string filename ) const
{
	//Read the POTCAR content
	std::string filecontent = this->get_textfile_constent( filename );

	const char * re = "(VRHFIN\\s*=\\s*\\w+\\s*:)";
	boost::regex expression(re);

	std::vector<std::string> blocks;
	std::copy(boost::sregex_token_iterator(filecontent.begin(), filecontent.end(), expression),
		boost::sregex_token_iterator(),
		std::back_inserter(blocks));

	const char * reAtom = "VRHFIN\\s*=\\s*(\\w+)\\s*:";
	boost::regex atom(reAtom);

	std::vector<std::string > result;
	for ( int ib = 0 ;ib < static_cast<int>(blocks.size()); ++ib)
	{
		boost::match_results<std::string::const_iterator> res;
		boost::regex_search(blocks[ib], res, atom );
		if ( res.size() != 2 )
			throw std::runtime_error(
					std::string("Failed to parse the atom list from POTCAR file ")+filename );
		result.push_back(std::string(res[1].first,res[1].second));
	}
	return result;
}

void VASPInterface::overwrite_POSCAR_file( std::string filename,
		std::vector<std::string > const & potcarAtomOrder,
		LatticeStructure::UnitCell const & unitcell ) const
{
	auto floatAccLine = [] (std::vector<double> v, int digits) {
		std::stringstream s;
		s << std::fixed << std::setprecision(digits);
		for ( auto x : v )
			s << x << " ";
		return s.str();
	};

	std::string fileContent;
	fileContent += "elephon created this file.\n";

	//scale factor
	fileContent += std::to_string(unitcell.get_alat())+"\n";

	//Lattice matrix
	auto A = unitcell.get_lattice_matrix();
	std::string a1str =
			floatAccLine( std::vector<double>( { A[3*0+0], A[3*1+0], A[3*2+0]} ) , 6 );
	fileContent += a1str+"\n";
	std::string a2str =
			floatAccLine( std::vector<double>( { A[3*0+1], A[3*1+1], A[3*2+1]} ) , 6 );
	fileContent += a2str+"\n";
	std::string a3str =
			floatAccLine( std::vector<double>( { A[3*0+2], A[3*1+2], A[3*2+2]} ) , 6 );
	fileContent += a3str+"\n";

	//atom types - here it gets tricky, because we need to match the order in the POTCAR file
	//The order is provided in potcarAtomOrder
	//First put the list of atom types ...
	for ( auto atom : potcarAtomOrder )
		fileContent += atom+" ";
	fileContent += "\n";

	// ... and how often they occur
	auto atoms = unitcell.get_atoms_list();
	std::map<std::string,int> numAtType;
	for ( auto a : atoms)
		if ( not numAtType.insert( std::make_pair(a.get_kind(),1) ).second )
			numAtType[ a.get_kind() ]++;

	for ( auto atom : potcarAtomOrder )
		fileContent += std::to_string( numAtType[atom] ) +" ";
	fileContent += "\n";

	//Now the coordinates
	fileContent += "direct\n";
	for ( auto atomFile : potcarAtomOrder )
		for ( auto a : atoms)
			if ( a.get_kind().compare( atomFile ) == 0 )//It is the turn of this atom kind
				fileContent += floatAccLine( a.get_position() , 6 )+" "+atomFile+"\n";

	std::ofstream file( filename.c_str() );
	file << fileContent;
}

void VASPInterface::write_KPOINTS_file(std::string filename,
		std::vector<double> const & kptShift,
		std::vector<int> const & monkhPackGrid) const
{
	assert( monkhPackGrid.size() == 3 );
	std::string fileContent = "elephon created this file\n"
			"0\n";
	bool gamma = true;
	for ( auto ks : kptShift )
		if ( gamma  and (std::abs(ks) < 1e-6 ))
			gamma = false;
	fileContent += gamma ? "Gamma Monkhorst-Pack\n" :  "Monkhorst-Pack\n";
	fileContent += std::to_string(monkhPackGrid[0])+" "
			+std::to_string(monkhPackGrid[1])+" "
			+std::to_string(monkhPackGrid[2])+"\n";
	fileContent += std::to_string(kptShift[0])+" "
			+std::to_string(kptShift[1])+" "
			+std::to_string(kptShift[2])+"\n";

	std::ofstream file( filename.c_str() );
	file << fileContent;
}

void VASPInterface::modify_incar_file(std::string filename,
		std::map<std::string,std::string> const & optionsToBeReset) const
{
	//Fortunately the scallop input file method can read the VASP INCAR file too ...
	IOMethods::InputFile inputfile;
	inputfile.read_input_file( filename );
	auto keyValuePairs = inputfile.get_all_input_config_values();

	//We have to insert the file options into the revised options, since the insert method will keep keys that match
	auto resultOptions = optionsToBeReset;
	resultOptions.insert( keyValuePairs.begin(), keyValuePairs.end() );

	std::string fileContent;
	for ( auto opt : resultOptions )
		fileContent += opt.first + " = " + opt.second + "\n";

	std::ofstream file( filename.c_str() );
	file << fileContent;
}

void VASPInterface::read_kpt_sampling(
		std::string root_directory,
		std::vector<int> & kptSampling,
		std::vector<double> & shifts)
{
	boost::filesystem::path root(root_directory);
	auto filename = (root / "KPOINTS").string();

	//make sure this is an automatic mesh
	std::ifstream file( filename.c_str() );
	std::string buffer;
	std::getline( file , buffer );//Comment
	std::getline( file , buffer );
	std::stringstream ss(buffer);
	int num = 1;
	ss >> num;
	if ( num != 0 )
		throw std::runtime_error( std::string("Can only use automatic k meshes in file ")+filename);

	kptSampling = std::vector<int>(3);
	shifts = std::vector<double>(3,0.0);

	std::getline( file , buffer );//Monkhorst-Pack flag

	std::getline( file , buffer );//Sampling
	std::stringstream ssKpts(buffer);
	ssKpts.exceptions( std::ios::failbit );
	for (int i = 0 ; i < 3 ; ++i)
		ssKpts >> kptSampling[i];

	std::getline( file , buffer );//shift - optional
	if ( buffer.empty() )
		return;
	std::stringstream ssKptshift(buffer);
	ssKptshift.exceptions( std::ios::failbit );
	for (int i = 0 ; i < 3 ; ++i)
		ssKpts >> shifts[i];

}

std::string VASPInterface::get_textfile_constent( std::string filename ) const
{
	std::ifstream file( filename.c_str() );
	file.seekg(0, std::ios::end);
	size_t size = file.tellg();
	std::string filecontent(size, ' ');
	file.seekg(0);
	file.read( &filecontent[0], size);
	return filecontent;
}

} /* namespace IOMethods */
} /* namespace elephon */
