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

#include "BuildFolderStructure.h"
#include <boost/filesystem.hpp>

namespace elephon
{
namespace IOMethods
{

void BuildFolderStructure::build(
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
	std::string sCellScfRoot = (elphd / "scell").string();
	std::vector<int> kptsSCell = input.get_kscell();
	for (int ki = 0 ; ki < static_cast<int>(kptsSCell.size()) ; ++ki )
		if ( kptsSCell[ki] <= 0)
			kptsSCell[ki] = std::ceil(double(defaultKSampling[ki])/double( input.get_scell()[ki] ));

	//Create the clean supercell
	elephon::LatticeStructure::UnitCell supercell =
			unitcell.build_supercell( input.get_scell().at(0), input.get_scell().at(1), input.get_scell().at(2));
	boost::filesystem::create_directories( elphd / "scell" );
	auto superCellOptions = interface.options_scf_supercell_no_wfctns_no_relax();

	interface.set_up_run(
			input.get_root_dir(),
			sCellScfRoot,
			kptsSCell,
			defaultKShift,
			supercell,
			superCellOptions);

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
