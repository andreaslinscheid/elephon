/*	This file main.cpp is part of elephon.
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
 *  Created on: 2017
 *      Author: A. Linscheid
 */
#include "IOMethods/Input.h"
#include "ElectronicStructure/LocalDensityOfStates.h"
#include "IOMethods/ChooseInterface.h"
#include "ElectronicStructure/BandStructureAnalysis.h"
#include "IOMethods/ResourceHandler.h"
#include "IOMethods/BuildFolderStructure.h"
#include "PhononStructure/AlphaSquaredF.h"

using namespace elephon;

int
main(int argc, char* argv[])
{
	IOMethods::Input input(argc,argv);
	auto options = input.get_opts();

	// Open the directory where we expect the preliminary calculation
	auto dataLoader = choose_interface(options);
	auto resources = std::make_shared<IOMethods::ResourceHandler>(dataLoader);

	// Print error messages if we are concerned the user is doing something stupid.
	dataLoader->check_prep_run( options.get_root_dir() );

	// Print the LDOS to file if desired ...
	if ( not options.get_f_ldos().empty() )
	{
		ElectronicStructure::LocalDensityOfStates ldos;
		ldos.compute_ldos(options.get_eldos(), resources);
		ldos.write_file(options.get_f_ldos());
	}

	// do auxiliary analysis of the band structure if desired ...
	ElectronicStructure::BandStructureAnalysis::do_band_structure_analysis(resources);

	// lep triggers the electron-phonon calculation
	if ( resources->get_optns().get_lep() or resources->get_optns().get_lp() )
	{
		IOMethods::BuildFolderStructure folderStructureBuilder;
		if ( folderStructureBuilder.check_is_build(resources->get_optns().get_elphd()) )
		{
			if ( resources->get_optns().get_lp() )
			{
				//run phonon calculation
				std::cout << "Calculating phonons ..." << std::endl;
				auto phGrid = resources->get_phonon_grid_obj();

				if ( not resources->get_optns().get_f_phdos().empty() )
				{
					std::cout << "Computing phonon DOS ..." << std::endl;
					// compute and write the phonon DOS
					auto frequencies = phGrid->setup_frequency_grid(
							resources->get_optns().get_phrange(),
							resources->get_optns().get_phnpts());
					phGrid->write_phonon_dos_file(
							resources->get_optns().get_f_phdos(),
							frequencies);
				}
			}

			if ( resources->get_optns().get_lep() )
			{
				//run electron-phonon calculation
				PhononStructure::AlphaSquaredF a2FDriver;
				std::cout << "Calculating a2F ..." << std::endl;
				a2FDriver.compute_a2F(resources);

				if ( not resources->get_optns().get_f_a2F().empty() )
					a2FDriver.write_a2F_file( resources->get_optns().get_f_a2F() );
			}
		}
		else
		{
			std::cout << "Building folder structure for phonon and electron-phonon calculation." << std::endl;
			folderStructureBuilder.build(resources);
		}
	}

    return 0;
}
