/*	This file AlphaSquaredF.cpp is part of elephon.
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
 *  Created on: Sep 26, 2017
 *      Author: A. Linscheid
 */

#include "PhononStructure/AlphaSquaredF.h"
#include "PhononStructure/ElectronPhononCoupling.h"
#include "ElectronicStructure/FermiSurface.h"
#include "Algorithms/FFTInterface.h"
#include <boost/filesystem.hpp>

namespace elephon
{
namespace PhononStructure
{

void
AlphaSquaredF::compute_a2F( std::shared_ptr<IOMethods::ElectronicStructureCodeInterface> dataLoader )
{
	std::string root_dir = dataLoader->get_optns().get_root_dir();
	boost::filesystem::path calcRoot( root_dir );
	std::string electrons_dir = (calcRoot / "electrons").string();

	// load the wavefunctions
	ElectronicStructure::Wavefunctions wfcts;
	wfcts.initialize( electrons_dir , dataLoader);

	// load the unit cell data
	LatticeStructure::UnitCell unitCell;
	dataLoader->read_unit_cell(electrons_dir, dataLoader->get_optns().get_gPrec(), unitCell);

	// load the band structure
	ElectronicStructure::ElectronicBands bands;
	dataLoader->read_band_structure(electrons_dir, bands);

	// rebuild the supercell.
	auto scd = dataLoader->get_optns().get_scell();
	auto supercell = unitCell.build_supercell(scd[0], scd[1], scd[2]);

	// rebuild the irreducible displacements
	std::vector<elephon::LatticeStructure::AtomDisplacement> irreducibleDispl;
	unitCell.generate_displacements(
			dataLoader->get_optns().get_magdispl(),
			dataLoader->get_optns().get_symDispl(),
			irreducibleDispl);

	// read the forces
	std::vector<std::vector<double>> forces(irreducibleDispl.size());
	for ( int irrep = 0 ; irrep < irreducibleDispl.size(); ++irrep)
	{
		auto irrepD = calcRoot / (std::string("displ_")+std::to_string(irrep));
		if ( ! boost::filesystem::exists(irrepD) )
			throw std::runtime_error(std::string("Error: Directory ")+irrepD.string()+" not present."
					" Must contain a converged VASP run");
		std::vector<double> f;
		dataLoader->read_forces(root_dir, f);
		forces[irrep]= std::move(f);
	}
	auto irrepPast = calcRoot / (std::string("displ_")+std::to_string(irreducibleDispl.size()));
	if ( boost::filesystem::exists(irrepPast) )
		throw std::runtime_error(std::string("Directory ")+irrepPast.string()+" is present but this run has only"
				+std::to_string(int(irreducibleDispl.size())-1)+" irreducible displacements.\n"
						"For safety reasons, the code is stopping here. Please clean up first.");

	// compute the phonon spectrum and displacement potential
	ForceConstantMatrix fc;
	fc.build(unitCell, supercell, irreducibleDispl, forces );
	Phonon ph;
	std::vector<double> masses;
	masses.reserve(unitCell.get_atoms_list().size());
	for ( auto const & a : unitCell.get_atoms_list() )
		masses.push_back(a.get_mass());
	ph.initialize(fc, masses);

	LatticeStructure::RegularBareGrid interpolationMesh;
	interpolationMesh.initialize(
			bands.get_grid().interpret_fft_dim_input(dataLoader->get_optns().get_fftd()),
			true,
			bands.get_grid().get_grid_prec(),
			dataLoader->get_optns().get_ffts(),
			bands.get_grid().get_lattice());

	// obtain a Fermi surface (set of constant energy surfaces) as a list of k point and weights
	int nSamples = dataLoader->get_optns().get_numFS();
	auto equalEnergySurfaces = dataLoader->get_optns().get_ea2f();
	for ( auto e : equalEnergySurfaces)
	{
		auto bndsCrossing = bands.get_bands_crossing_energy_lvls({e});
		std::vector<double> reducibleData;
		bands.generate_interpolated_reducible_grid_bands(
				bndsCrossing,
				interpolationMesh,
				reducibleData);

		ElectronicStructure::FermiSurface fs;
		fs.triangulate_surface(
				interpolationMesh,
				bndsCrossing.size(),
				reducibleData,
				nSamples,
				e);


	// integrate the function on each constant energy surfaces.

	ElectronPhononCoupling gkkp;
//	gkkp.generate_gkkp( kpts, kppts, bandsList, bandspList, ph, dvscf, wfcts );
	}
}

} /* namespace PhononStructure */
} /* namespace elephon */
