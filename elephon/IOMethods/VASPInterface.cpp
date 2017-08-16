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

std::string
VASPInterface::code_tag() const
{
	return "vasp";
}

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

	//first things first, we copy the POTCAR file to the location
	boost::filesystem::path potcarPrev = root / "POTCAR" ;
	boost::filesystem::path potcarNew = elphd /  "POTCAR";
	boost::filesystem::copy( potcarPrev, potcarNew );

	//write parameters in POSCAR according to data in unitcell
	auto atomOrder = this->read_potcar_atom_order( potcarNew.string() );
	boost::filesystem::path poscarNew = elphd / "POSCAR";
	this->overwrite_POSCAR_file( poscarNew.string(), atomOrder, unitcell );

	//Write the KPOINTS file
	boost::filesystem::path kpts = elphd / "KPOINTS";
	this->write_KPOINTS_file( kpts.string(), kptShift, kptSampling );

	//First copy and then modify the INCAR file
	boost::filesystem::path incarPrev = root / "INCAR" ;
	boost::filesystem::path incarNew = elphd /  "INCAR";
	boost::filesystem::copy( incarPrev, incarNew );

	this->modify_incar_file( incarNew.string(), options );
}

std::map<std::string,std::string>
VASPInterface::options_nscf_keep_wfctns_no_relax() const
{
	std::map<std::string,std::string> options;
	//We want to keep the wavefunctions
	options["LWAVE"] = ".TRUE.";
	options["IBRION"] = "-1";
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
	options["IBRION"] = "-1";
	options["NSW"] = "0";
	options["ICHARG"] = "1";
	options["LVTOT"] = ".TRUE.";
	options["PREC"] = "High";
	options["LVHAR"] = ".FALSE.";
	return options;
}

void
VASPInterface::read_wavefunctions(
		std::string root_directory,
		std::vector<int> const & kpts,
		std::vector<int> const & bandIndices,
		std::vector< std::vector< std::complex<float> > > & wfctData,
		std::vector< int > & npwPerKpt)
{
	boost::filesystem::path rootdir(root_directory);
	wfcReader_.prepare_wavecar( (rootdir / "WAVECAR").string() );
	wfcReader_.read_wavefunction( kpts, bandIndices, wfctData, npwPerKpt );
}

void
VASPInterface::compute_fourier_map(
		std::vector<double> const & kpts,
		std::vector< std::vector<int> > & fourierMap,
		double gridPrec)
{
	if ( wfcReader_.get_filename().empty() )
		throw std::runtime_error("Cannot compute Fourier map without pointing to a WAVECAR file that specifies the cutoff");
	//VASP does not store these mapping on disk.
	wfcReader_.compute_fourier_map(kpts, fourierMap, gridPrec);
}

std::vector<int>
VASPInterface::get_max_fft_dims() const
{
	return wfcReader_.get_fft_max_dims();
}

void
VASPInterface::read_cell_paramters(
		std::string root_directory,
		double symPrec,
		LatticeStructure::RegularSymmetricGrid & kPointMesh,
		LatticeStructure::LatticeModule & lattice,
		std::vector<LatticeStructure::Atom> & atoms,
		LatticeStructure::Symmetry & symmetry)
{
	this->check_open_poscar(root_directory);
	lattice.initialize( posReader_.get_lattice_matrix() );
	atoms = posReader_.get_atoms_list();
	if ( symReader_.get_symmetries().empty() )
		symReader_.read_file( (boost::filesystem::path(root_directory) / "OUTCAR").string() );
	symmetry.initialize( symPrec, symReader_.get_symmetries(),
			symReader_.get_fractionTranslations(), lattice, symReader_.get_time_revesal_symmetry() );

	//Read the k point in the irreducible zone
	std::vector<int> kDim;
	std::vector<double> shifts;
	this->read_kpt_sampling(root_directory, kDim, shifts);
	LatticeStructure::Symmetry kpointSymmetry = symmetry;

	//make sure that the irreducible zone of elephon and VASP agree
	boost::filesystem::path rootdir(root_directory);
	std::vector<double> irreducibleKPoints;
	if ( boost::filesystem::exists( rootdir / "WAVECAR" ) )
	{
		wfcReader_.prepare_wavecar( (rootdir / "WAVECAR").string() );
		irreducibleKPoints = wfcReader_.get_k_points();
	}
	else
	{
		xmlReader_.parse_file(  (rootdir / "vasprun.xml").string() );
		irreducibleKPoints = xmlReader_.get_k_points();
	}

	kpointSymmetry.set_reciprocal_space_sym();
	kPointMesh.initialize( kDim, symPrec, shifts, kpointSymmetry, lattice, irreducibleKPoints );
}

void
VASPInterface::read_lattice_structure(
		std::string root_directory,
		LatticeStructure::LatticeModule & lattice)
{
	this->check_open_poscar(root_directory);
	lattice.initialize( posReader_.get_lattice_matrix() );
}

void
VASPInterface::check_open_poscar(std::string const & root_dir )
{
	boost::filesystem::path rootdir(root_dir);
	//See if we find a POTCAR file which supplies the atom names
	auto atomOrder = read_potcar_atom_order((rootdir / "POTCAR").string());
	if ( posReader_.get_atoms_list().empty() )
		posReader_.read_file( (rootdir / "POSCAR").string(), atomOrder );
}

void
VASPInterface::read_forces(
		std::string root_directory,
		std::vector<double> & forces)
{
	boost::filesystem::path rootdir(root_directory);
	xmlReader_.parse_file( (rootdir / "vasprun.xml").string() );
	forces = xmlReader_.get_forces();
}

void VASPInterface::read_electronic_potential(
		std::string root_directory,
		std::vector<int> & dims,
		std::vector<double> & output)
{
	boost::filesystem::path rootdir(root_directory);
	potReader_.read_scf_potential( (rootdir / "LOCPOT").string(), dims, output );
}

void
VASPInterface::read_band_structure(
		std::string root_directory,
		ElectronicStructure::ElectronicBands & bands)
{
	//Load the irreducible data
	int nkIrred = 0;
	int nband = 0;
	double eFermi = 0.0;
	std::vector<double> irredData;
	this->read_electronic_structure(
			root_directory,
			nband,
			nkIrred,
			irredData,
			eFermi);

	//expand to the reducible zone
	LatticeStructure::RegularSymmetricGrid kPointMesh;
	LatticeStructure::LatticeModule lattice;
	std::vector<LatticeStructure::Atom> atoms;
	LatticeStructure::Symmetry symmetry;
	this->read_cell_paramters(
			root_directory,
			this->get_optns().get_gPrec(),
			kPointMesh,
			lattice,
			atoms,
			symmetry);

	if ( kPointMesh.get_np_irred() != nkIrred )
		throw std::runtime_error("Number of k point in the energy grid and the loaded k point grid disagree.");

	bands.initialize( nband, irredData, kPointMesh );
}

void
VASPInterface::read_electronic_structure(
		std::string root_directory,
		int & nBnd,
		int & nkptsIrred,
		std::vector<double> & energies,
		double & fermiEnergy)
{
	boost::filesystem::path rootdir(root_directory);
	xmlReader_.parse_file( (rootdir / "vasprun.xml").string() );
	fermiEnergy = xmlReader_.get_Fermi_energy();

	wfcReader_.prepare_wavecar( (rootdir / "WAVECAR").string() );
	energies = wfcReader_.get_energies();

	nBnd = wfcReader_.get_num_bands();
	nkptsIrred = wfcReader_.get_num_kpts();
}


std::vector<std::string >
VASPInterface::read_potcar_atom_order( std::string filename ) const
{
	//Read the POTCAR content
	std::string filecontent = this->get_textfile_content( filename );

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
	std::string a1str =
			floatAccLine( unitcell.get_lattice().get_lattice_vector(0) , 6 );
	fileContent += a1str+"\n";
	std::string a2str =
			floatAccLine( unitcell.get_lattice().get_lattice_vector(1) , 6 );
	fileContent += a2str+"\n";
	std::string a3str =
			floatAccLine( unitcell.get_lattice().get_lattice_vector(2) , 6 );
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
			{
				//map to the zone [0,1[
				auto pos = a.get_position();
				for ( auto & xi : pos )
					xi = xi < 0 ? xi + 1.0 : xi;

				fileContent += floatAccLine( pos , 6 )+" "+atomFile+"\n";
			}

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
		std::vector<double> & shifts) const
{
	boost::filesystem::path root(root_directory);
	auto filename = (root / "KPOINTS").string();

	//make sure this is an automatic mesh
	std::ifstream file( filename.c_str() );
	if ( ! file.good() )
		throw std::runtime_error( filename + ": file not readable" );
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
	bool isGammaSet = false;
	if ( ! buffer.empty() )
		isGammaSet = ((buffer.front()  == 'G') || (buffer.front()  == 'g'));

	std::getline( file , buffer );//Sampling
	std::stringstream ssKpts(buffer);
	ssKpts.exceptions( std::ios::failbit );
	for (int i = 0 ; i < 3 ; ++i)
		ssKpts >> kptSampling[i];

	std::getline( file , buffer );//shift - optional
	if ( buffer.empty() )
	{
		//For even meshes shift by half a cell
		for ( int  i = 0 ; i < 3 ; ++i )
			if ( (kptSampling[i]%2 == 0) and (not isGammaSet) )
				shifts[i] = 0.5;
		return;
	}
	std::stringstream ssKptshift(buffer);
	ssKptshift.exceptions( std::ios::failbit );
	for (int i = 0 ; i < 3 ; ++i)
		ssKptshift >> shifts[i];

}

std::string VASPInterface::get_textfile_content( std::string filename ) const
{
	std::ifstream file( filename.c_str() );
	if ( not file.good() )
		throw std::runtime_error( std::string("file ")+filename+ " not readable" );
	file.seekg(0, std::ios::end);
	size_t size = file.tellg();
	std::string filecontent(size, ' ');
	file.seekg(0);
	file.read( &filecontent[0], size);
	return filecontent;
}

} /* namespace IOMethods */
} /* namespace elephon */
