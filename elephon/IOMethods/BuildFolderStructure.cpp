/*	This file BuildFolderStructure.cpp is part of elephon.
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

#include "IOMethods/BuildFolderStructure.h"
#include <boost/filesystem.hpp>
#include <stdexcept>

namespace elephon
{
namespace IOMethods
{

void
BuildFolderStructure::build(std::shared_ptr<ResourceHandler> resources) const
{
	boost::filesystem::path elphd(resources->get_optns().get_elphd());
	if ( this->check_is_build(elphd.string()) )
		throw std::runtime_error( std::string("Refusing to overwrite ")
				+elphd.string()+" - delete it to regenerate." );

	auto interface = resources->get_electronic_structure_interface();
	auto kgrid = resources->get_electronic_bands_obj()->get_grid();
	auto unitcell = resources->get_primitive_unitcell_obj();

	if ( not resources->get_optns().get_eld().empty() )
	{
		//Create the run files for the dense electronic structure
		std::string denseNscfRoot = resources->get_optns().get_eld();
		boost::filesystem::create_directories( denseNscfRoot );

		auto kptsDense = kgrid.interpret_fft_dim_input( resources->get_optns().get_kdense() );
		auto kShift = resources->get_optns().get_ksdense();

		auto denseNSCFOptions = interface->options_nscf_keep_wfctns_no_relax();
		interface->set_up_run(
				resources->get_optns().get_root_dir(),
				denseNscfRoot,
				kptsDense,
				kShift,
				*unitcell,
				denseNSCFOptions);

		interface->copy_charge(
				resources->get_optns().get_root_dir(),
				denseNscfRoot);
	}

	//Determine the kpt sampling for the supercells
	std::vector<int> kptsSCell = resources->get_optns().get_kscell();
	for (int ki = 0 ; ki < static_cast<int>(kptsSCell.size()) ; ++ki )
		if ( kptsSCell[ki] <= 0)
			kptsSCell[ki] = std::ceil(double(kgrid.get_grid_dim()[ki])
									 /double(resources->get_optns().get_scell()[ki]) );

	//Create the clean supercell
	std::vector<int> scMulti{	resources->get_optns().get_scell().at(0),
								resources->get_optns().get_scell().at(1),
								resources->get_optns().get_scell().at(2)};
	elephon::LatticeStructure::UnitCell supercell =
			unitcell->build_supercell(scMulti[0], scMulti[1], scMulti[2]);
	if ( supercell.get_symmetry().get_num_symmetries() != unitcell->get_symmetry().get_num_symmetries() )
		throw std::runtime_error("Unable to handle symmetry reduction duing supercell create : \n"
				"make the supercell respect the unit cell symmetry!");

	auto superCellOptions = interface->options_scf_supercell_no_wfctns_no_relax();
	// set the real space grid dimensions explicitly
	auto rswfctDim = interface->read_wfct_real_space_grid_dim(resources->get_optns().get_root_dir());
	superCellOptions["NGX"] = std::to_string(rswfctDim[0] * scMulti[0]);
	superCellOptions["NGY"] = std::to_string(rswfctDim[1] * scMulti[1]);
	superCellOptions["NGZ"] = std::to_string(rswfctDim[2] * scMulti[2]);
	auto rsChargeDim = interface->read_charge_real_space_grid_dim(resources->get_optns().get_root_dir());
	superCellOptions["NGXF"] = std::to_string(rsChargeDim[0] * scMulti[0]);
	superCellOptions["NGYF"] = std::to_string(rsChargeDim[1] * scMulti[1]);
	superCellOptions["NGZF"] = std::to_string(rsChargeDim[2] * scMulti[2]);

	// create displacements of atoms in the primitive cell and apply them in the supercell
	std::vector< elephon::LatticeStructure::AtomDisplacement>  irreducible;
	unitcell->generate_displacements(
			resources->get_optns().get_magdispl(),
			resources->get_optns().get_symDispl(),
			irreducible );

	//transform them from the primitive to the supercell
	double s[3] = {1.0/scMulti[0], 1.0/scMulti[1], 1.0/scMulti[2]};
	for ( auto &id : irreducible )
		id.scale_position(s[0], s[1], s[2]);

	//We only care about the irreducible displacement at this point.
	for ( int ipert = 0 ; ipert < static_cast<int>(irreducible.size()); ipert++ )
	{
		auto perti = supercell;
		perti.displace_atom( irreducible[ipert] );
		boost::filesystem::path dir = elphd / (std::string("displ_")+std::to_string(ipert));
		boost::filesystem::create_directories( dir );

		interface->set_up_run(
				resources->get_optns().get_root_dir(),
				dir.string(),
				kptsSCell,
				kgrid.get_grid_shift(),
				perti,
				superCellOptions);
	}
}

void
BuildFolderStructure::build(
		IOMethods::InputOptions const & input,
		LatticeStructure::UnitCell const & unitcell,
		ElectronicStructureCodeInterface & interface ) const
{
	boost::filesystem::path elphd( input.get_elphd() );
	if ( boost::filesystem::exists( elphd ) )
		throw std::runtime_error( std::string("Refusing to overwrite "
				+elphd.string()+" - delete it to regenerate."));

	//Set up the wavefunction/energies grid - see if we determine the sampling automatically
	const int scale_factor = 4;
	std::string denseNscfRoot = (elphd / "electrons").string();
	std::vector<int> defaultKSampling;
	std::vector<double> defaultKShift;
	interface.read_kpt_sampling( input.get_root_dir(),  defaultKSampling, defaultKShift );
	std::vector<int> kptsDense = input.get_kdense();
	for (int ki = 0 ; ki < static_cast<int>(kptsDense.size()) ; ++ki )
		if ( kptsDense[ki] <= 0)
			kptsDense[ki] = defaultKSampling[ki]*scale_factor;

	//Create the run files for the electronic structure
	auto denseNSCFOptions = interface.options_nscf_keep_wfctns_no_relax();
	boost::filesystem::create_directories( denseNscfRoot );
	interface.set_up_run(
			input.get_root_dir(),
			denseNscfRoot,
			kptsDense,
			defaultKShift,
			unitcell,
			denseNSCFOptions);

	//Determine the kpt sampling for the supercells
	std::vector<int> kptsSCell = input.get_kscell();
	for (int ki = 0 ; ki < static_cast<int>(kptsSCell.size()) ; ++ki )
		if ( kptsSCell[ki] <= 0)
			kptsSCell[ki] = std::ceil(double(defaultKSampling[ki])/double( input.get_scell()[ki] ));
	auto superCellOptions = interface.options_scf_supercell_no_wfctns_no_relax();

	//Create the clean supercell
	elephon::LatticeStructure::UnitCell supercell =
			unitcell.build_supercell( input.get_scell().at(0), input.get_scell().at(1), input.get_scell().at(2));

	//Create displacements - we have to reduce the symmetry to the supercell in case this is different
	LatticeStructure::UnitCell unitcellReducedSym;
	unitcellReducedSym.initialize( unitcell.get_atoms_list(), unitcell.get_lattice(), supercell.get_symmetry() );
	std::vector< elephon::LatticeStructure::AtomDisplacement>  irreducible;
	unitcellReducedSym.generate_displacements(
			input.get_magdispl(),
			input.get_symDispl(),
			irreducible );

	//transform them from the primitive to the supercell
	double s[3] = {1.0/input.get_scell().at(0) , 1.0/input.get_scell().at(1), 1.0/input.get_scell().at(2)};
	for ( auto &id : irreducible )
		id.scale_position( s[0], s[1], s[2] );

	//We only care about the irreducible displacement at this point.
	for ( int ipert = 0 ; ipert < static_cast<int>(irreducible.size()); ipert++ )
	{
		auto perti = supercell;
		perti.displace_atom( irreducible[ipert] );
		boost::filesystem::path dir = elphd / (std::string("displ_")+std::to_string(ipert));
		boost::filesystem::create_directories( dir );

		interface.set_up_run(
				input.get_root_dir(),
				dir.string(),
				kptsSCell,
				defaultKShift,
				perti,
				superCellOptions);
	}
}

bool BuildFolderStructure::check_is_build( std::string rootFolder ) const
{
	boost::filesystem::path root_dir(rootFolder);
	if ( not boost::filesystem::exists(root_dir ) )
		return false;
	return true;
}

} /* namespace IOMethods */
} /* namespace elephon */
