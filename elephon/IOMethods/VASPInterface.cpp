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
#include "IOMethods/ReadVASPxmlFile.h"
#include "AtomicSite/AtomSiteData.h"
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <map>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <algorithm>

namespace elephon
{
namespace IOMethods
{

std::string
VASPInterface::code_tag() const
{
	return "vasp";
}

void
VASPInterface::check_prep_run(
		std::string root_directory ) const
{
	//Fortunately the scallop input file method can read the VASP INCAR file too ...
	boost::filesystem::path root(root_directory);
	IOMethods::InputFile inputfile;
	inputfile.read_input_file( (root / "INCAR").string() );
	auto keyValuePairs = inputfile.get_all_input_config_values();

	auto print_warning = [] (std::string const & message) {
		std::string completeMessage = std::string("WARNING: ")+message+"\n";
		int nEqalChars = 0;
		std::stringstream ss(completeMessage);
		std::string item;
		while (std::getline(ss, item, '\n'))
		{
			nEqalChars = std::max(nEqalChars, static_cast<int>(item.size()));
		}
		std::cerr << std::string(nEqalChars, '=') << "\n"
				  << completeMessage;
		std::cerr << std::string(nEqalChars, '=') <<std::endl;
	};

	// Check if the user is running the tetrahedra integration scheme which is not a good idea
	// when phonons / electron-phonon calculations are done.
	if ( this->get_optns().get_lep() or this->get_optns().get_lp())
	{
		auto ret = keyValuePairs.find("ISMEAR");
		if ( ret != keyValuePairs.end() )
			if ( std::stoi(ret->second) < 0 )
				print_warning("It seems you are using the Tetrahedra BZ integration scheme for a phonon calculation\n"
						"This is not a good idea (see VASP manual). Continue only if you know what you are doing!!");
	}

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
	std::vector<std::string> atomNames(atomOrder.size());
	for ( int i = 0 ; i < atomOrder.size(); ++i)
		atomNames[i] = atomOrder[i].first;
	boost::filesystem::path poscarNew = elphd / "POSCAR";
	this->overwrite_POSCAR_file( poscarNew.string(), atomNames, unitcell );

	//Write the KPOINTS file
	boost::filesystem::path kpts = elphd / "KPOINTS";
	this->write_KPOINTS_file( kpts.string(), kptShift, kptSampling );

	//First copy and then modify the INCAR file
	boost::filesystem::path incarPrev = root / "INCAR" ;
	boost::filesystem::path incarNew = elphd /  "INCAR";
	boost::filesystem::copy( incarPrev, incarNew );

	this->modify_incar_file( incarNew.string(), options );
}

void
VASPInterface::copy_charge(
		std::string root_directory,
		std::string target_directory) const
{
	boost::filesystem::path root(root_directory);
	boost::filesystem::path elphd(target_directory);
	boost::filesystem::path chgcarPrev = root / "CHGCAR" ;
	boost::filesystem::path chgcarNew = elphd /  "CHGCAR";
	boost::filesystem::copy( chgcarPrev, chgcarNew );
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
	options["LVHAR"] = ".FALSE.";
	options["LCHARG"] = ".FALSE.";
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
VASPInterface::get_max_fft_dims()
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
	boost::filesystem::path rootdir(root_directory);
	this->read_lattice_structure(root_directory, lattice);
	this->read_atoms_list(root_directory, atoms);
	this->read_symmetry(root_directory, symPrec, lattice, symmetry);

	//Read the k point in the irreducible zone
	std::vector<int> kDim;
	std::vector<double> shifts;
	this->read_kpt_sampling(root_directory, kDim, shifts);
	LatticeStructure::Symmetry kpointSymmetry = symmetry;

	//make sure that the irreducible zone of elephon and VASP agree
	std::vector<double> irreducibleKPoints;
	if ( boost::filesystem::exists( rootdir / "WAVECAR" ) )
	{
		wfcReader_.prepare_wavecar( (rootdir / "WAVECAR").string() );
		irreducibleKPoints = wfcReader_.get_k_points();
	}
	else if ( boost::filesystem::exists( rootdir / "vasprun.xml" ) )
	{
		xmlReader_.parse_file(  (rootdir / "vasprun.xml").string() );
		irreducibleKPoints = xmlReader_.get_k_points();
	}

	kpointSymmetry.set_reciprocal_space_sym();
	kPointMesh.initialize( kDim, symPrec, shifts, kpointSymmetry, lattice, irreducibleKPoints );
}

std::vector<int>
VASPInterface::read_wfct_real_space_grid_dim(std::string root_directory)
{
	boost::filesystem::path rootdir(root_directory);
	xmlReader_.parse_file(  (rootdir / "vasprun.xml").string() );
	return xmlReader_.get_wfct_fourier_dim();
}

std::vector<int>
VASPInterface::read_charge_real_space_grid_dim(std::string root_directory)
{
	boost::filesystem::path rootdir(root_directory);
	xmlReader_.parse_file(  (rootdir / "vasprun.xml").string() );
	return xmlReader_.get_charge_fourier_dim();
}


void
VASPInterface::read_unit_cell(
		std::string root_directory,
		double symprec,
		LatticeStructure::UnitCell & unitcell )
{
	boost::filesystem::path rootdir(root_directory);
	LatticeStructure::LatticeModule lattice;
	this->read_lattice_structure(root_directory, lattice);
	std::vector<LatticeStructure::Atom> atoms;
	this->read_atoms_list(root_directory, atoms);
	LatticeStructure::Symmetry symmetry;
	this->read_symmetry(root_directory, symprec, lattice, symmetry);
	unitcell.initialize(atoms, lattice, symmetry);
}

void
VASPInterface::read_lattice_structure(
		std::string root_directory,
		LatticeStructure::LatticeModule & lattice)
{
	boost::filesystem::path rootdir(root_directory);
	if ( not boost::filesystem::exists(rootdir) )
		throw std::runtime_error(std::string("Directory ")+root_directory + " does not exist."
				" Failed to read lattice structure.");

	if (boost::filesystem::exists(rootdir / "vasprun.xml" ))
	{
		xmlReader_.parse_file( (rootdir / "vasprun.xml").string() );
		lattice.initialize( xmlReader_.get_lattice_matrix() );
		return;
	}

	if ( boost::filesystem::exists(rootdir / "POSCAR" ) )
	{
		if ( posReader_.get_atoms_list().empty() )
		{
			if ( ! boost::filesystem::exists(rootdir / "POTCAR" ) )
				throw std::runtime_error("Cannot complete atom data without data from POTCAR");
			auto atomOrder = read_potcar_atom_order((rootdir / "POTCAR").string());
			posReader_.read_file( (rootdir / "POSCAR").string(), atomOrder );
		}
		lattice.initialize( posReader_.get_lattice_matrix() );
		return;
	}

	throw std::runtime_error(std::string("No file to parse structure from.\n") +
				" Need either a POSCAR or a vasprun.xml file in directory "+root_directory);
}

void
VASPInterface::read_atoms_list(
		std::string root_directory,
		std::vector<LatticeStructure::Atom> & atoms)
{
	boost::filesystem::path rootdir(root_directory);

	if (boost::filesystem::exists(rootdir / "vasprun.xml" ))
	{
		xmlReader_.parse_file( (rootdir / "vasprun.xml").string() );
		atoms = xmlReader_.get_atoms_list();
		return;
	}

	if ( boost::filesystem::exists(rootdir / "POSCAR" ) )
	{
		if ( posReader_.get_atoms_list().empty() )
		{
			if ( ! boost::filesystem::exists(rootdir / "POTCAR" ) )
				throw std::runtime_error("Cannot complete atom data without data from POTCAR");
			auto atomOrder = read_potcar_atom_order((rootdir / "POTCAR").string());
			posReader_.read_file( (rootdir / "POSCAR").string(), atomOrder );
		}
		atoms = posReader_.get_atoms_list();
		return;
	}

	throw std::runtime_error("No file to parse structure from");
}

void
VASPInterface::read_symmetry(
		std::string root_directory,
		double symprec,
		LatticeStructure::LatticeModule const& lattice,
		LatticeStructure::Symmetry & symmetry)
{
	boost::filesystem::path rootdir(root_directory);

	if ( symReader_.get_symmetries().empty() )
		symReader_.read_file( (rootdir / "OUTCAR").string() );
	symmetry.initialize(
			symprec,
			symReader_.get_symmetries(),
			symReader_.get_fractionTranslations(),
			lattice,
			symReader_.get_time_revesal_symmetry() );
}

void
VASPInterface::check_open_poscar(std::string const & root_dir )
{
	boost::filesystem::path rootdir(root_dir);
	//See if we find a POTCAR file which supplies the atom names
	if ( posReader_.get_atoms_list().empty() )
	{
		auto atomOrder = read_potcar_atom_order((rootdir / "POTCAR").string());
		posReader_.read_file( (rootdir / "POSCAR").string(), atomOrder );
	}
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
		Auxillary::alignedvector::DV & regularGrid,
		std::vector<AtomicSite::AtomSiteData> & radialPart)
{
	boost::filesystem::path rootdir(root_directory);
	std::vector<double> regularGridPartV;
	potReader_.read_scf_potential( (rootdir / "LOCPOT").string(), dims, regularGridPartV );
	Auxillary::alignedvector::DV regularGridPart(regularGridPartV.begin(), regularGridPartV.end());
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

	bands.initialize(
			nband,
			eFermi,
			std::move(Auxillary::alignedvector::DV(irredData.begin(), irredData.end())),
			kPointMesh );
}

void
VASPInterface::read_nBnd(
		std::string root_directory,
		int & nBnd)
{
	boost::filesystem::path rootdir(root_directory);
	// this information is either in the wavecar or in the vasprun.xml
	if ( boost::filesystem::exists( rootdir / "WAVECAR" ) )
	{
		wfcReader_.prepare_wavecar( (rootdir / "WAVECAR").string() );
		nBnd = wfcReader_.get_num_bands();
	}
	else if ( boost::filesystem::exists( rootdir / "vasprun.xml" ) )
	{
		xmlReader_.parse_file( (rootdir / "vasprun.xml").string() );
		nBnd = xmlReader_.get_nBnd();
	}
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

	energies = xmlReader_.get_energies();

	nBnd = xmlReader_.get_nBnd();
	nkptsIrred = xmlReader_.get_nkp();
}


std::vector<std::pair<std::string, double> >
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

	std::vector<std::string > atoms;
	for ( int ib = 0 ;ib < static_cast<int>(blocks.size()); ++ib)
	{
		boost::match_results<std::string::const_iterator> res;
		boost::regex_search(blocks[ib], res, atom );
		if ( res.size() != 2 )
			throw std::runtime_error(
					std::string("Failed to parse the atom list from POTCAR file ")+filename );
		atoms.push_back(std::string(res[1].first,res[1].second));
	}

	std::vector<std::pair<std::string, double> > result(atoms.size());

	const char * re_mass = "(POMASS\\s*=\\s*\\d+\\.\\d+\\s*;)";
	boost::regex expression_mass(re_mass);

	std::vector<std::string> blocks_mass;
	std::copy(boost::sregex_token_iterator(filecontent.begin(), filecontent.end(), expression_mass),
		boost::sregex_token_iterator(),
		std::back_inserter(blocks_mass));

	if ( (blocks_mass.size() != atoms.size()) or (blocks_mass.size() == 0) )
		throw std::runtime_error(std::string("Error parsing POTCAR ")+filename+ " for atom names and mass");

	const char * reAtomMass = "POMASS\\s*=\\s*(\\d+\\.\\d+)\\s*;";
	boost::regex mass(reAtomMass);
	for ( int ib = 0 ;ib < static_cast<int>(blocks_mass.size()); ++ib)
	{
		boost::match_results<std::string::const_iterator> res;
		boost::regex_search(blocks_mass[ib], res, mass );
		if ( ib >= result.size() )
			throw std::runtime_error(
					std::string("Failed to parse the atom mass; incompatible mass and names in POTCAR file ")+filename );

		if ( res.size() != 2 )
			throw std::runtime_error(
					std::string("Failed to parse the atom mass list from POTCAR file ")+filename );
		result[ib] = std::move(std::make_pair(atoms[ib], std::stof(std::string(res[1].first,res[1].second)) ));
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
	fileContent += " 1.0\n";

	auto a1 = unitcell.get_lattice().get_lattice_vector(0);
	for (auto & xi : a1 )
		xi *= unitcell.get_alat();

	auto a2 = unitcell.get_lattice().get_lattice_vector(1);
	for (auto & xi : a2 )
		xi *= unitcell.get_alat();

	auto a3 = unitcell.get_lattice().get_lattice_vector(2);
	for (auto & xi : a3 )
		xi *= unitcell.get_alat();

	//Lattice matrix
	std::string a1str =	floatAccLine( a1 , 12 );
	fileContent += a1str+"\n";
	std::string a2str = floatAccLine( a2 , 12 );
	fileContent += a2str+"\n";
	std::string a3str = floatAccLine( a3 , 12 );
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

				fileContent += floatAccLine( pos , 12 )+" "+atomFile+"\n";
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
	fileContent += "Gamma Monkhorst-Pack\n"; // internally, the k point shift is always measured from zero.
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
	if ( boost::filesystem::exists( root / "vasprun.xml" ) )
	{
		xmlReader_.parse_file( (root / "vasprun.xml").string() );
		kptSampling = xmlReader_.get_k_grid_dim();
		shifts = xmlReader_.get_k_grid_shift();
		return;
	}

	if ( boost::filesystem::exists( (root / "KPOINTS") ) )
	{

		LatticeStructure::LatticeModule lattice;
		this->read_lattice_structure(root.string(), lattice);
		kpointReader_.read_kpoints((root / "KPOINTS").string(), lattice);

		kptSampling = kpointReader_.get_grid_dim();
		shifts = kpointReader_.get_grid_shift();
		return;
	}

	throw std::runtime_error("Need vasprun.xml or KPOINTS to read the k grid parameters");
}

void
VASPInterface::read_reciprocal_symmetric_grid(
		std::string root_directory,
		LatticeStructure::RegularSymmetricGrid & kgrid)
{
	boost::filesystem::path rootdir(root_directory);
	std::vector<int> kptSampling;
	std::vector<double> shifts;
	this->read_kpt_sampling(root_directory, kptSampling, shifts);

	LatticeStructure::LatticeModule lattice;
	this->read_lattice_structure(root_directory, lattice);

	LatticeStructure::Symmetry sym;
	this->read_symmetry(root_directory, this->get_optns().get_gPrec(), lattice, sym);
	sym.set_reciprocal_space_sym(true);

	//make sure that the irreducible zone of elephon and VASP agree
	std::vector<double> irreducibleKPoints;
	if ( boost::filesystem::exists( rootdir / "WAVECAR" ) )
	{
		wfcReader_.prepare_wavecar( (rootdir / "WAVECAR").string() );
		irreducibleKPoints = wfcReader_.get_k_points();
	}
	else if ( boost::filesystem::exists( rootdir / "vasprun.xml" ) )
	{
		xmlReader_.parse_file(  (rootdir / "vasprun.xml").string() );
		irreducibleKPoints = xmlReader_.get_k_points();
	}
	LatticeStructure::RegularSymmetricGrid g;
	g.initialize(kptSampling, this->get_optns().get_gPrec(), shifts, sym, lattice, irreducibleKPoints);
	kgrid = g;
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
