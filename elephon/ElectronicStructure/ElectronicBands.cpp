/*	This file ElectronicBands.cpp is part of elephon.
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
 *  Created on: May 20, 2017
 *      Author: A. Linscheid
 */

#include "ElectronicBands.h"
#include "Algorithms/FFTInterface.h"
#include "Algorithms/LinearAlgebraInterface.h"
#include <stdexcept>
#include <set>
#include <chrono>
#include <fstream>

namespace elephon
{
namespace ElectronicStructure
{

void
ElectronicBands::initialize(
		std::vector<double> const & kpoints,
		int numBands,
		double fermiEnergy,
		std::vector<double> bandData,
		LatticeStructure::RegularSymmetricGrid grid)
{
	assert( grid.is_reci() );
	std::vector<int> pointIndices;
	grid.get_list_lattice_point_indices(kpoints, pointIndices);

	//Check for errors
	int npIrr = pointIndices.size();
	if( not ( npIrr ==  grid.get_np_irred() ) )
		throw std::runtime_error("Initialization data must be complete in the irreducible zone.");

	if ( bandData.size() != static_cast<size_t>(numBands)*kpoints.size()/3 )
		throw std::logic_error("Called initialize inconsistent data size.");

	std::vector<double> dataIrredOrdered(npIrr*numBands);
	for ( int irr = 0 ; irr < npIrr; ++irr )
		for ( int ibnd = 0 ; ibnd < numBands; ++ibnd )
			dataIrredOrdered[ pointIndices[irr]*numBands + ibnd ] = bandData[ irr*numBands + ibnd ] - fermiEnergy;

	LatticeStructure::DataRegularGrid<double>::initialize(
			numBands,
			fermiEnergy,
			std::move(dataIrredOrdered),
			std::move(grid));
}

int
ElectronicBands::get_nspin() const
{
	return 1;
}

int
ElectronicBands::get_nBnd() const
{
	return this->get_nData_gpt()/this->get_nspin();
}

double
ElectronicBands::operator() (int ikIrred, int ib, int ispin) const
{
	assert( (ispin >= 0) && (ispin < this->get_nspin()));
	return read(ikIrred, ib*this->get_nspin()+ispin);
}

ElectronicBands
ElectronicBands::fft_interpolate_part(
		int startBnD, int endBnd,
		std::vector<int> const & newDims,
		std::vector<double> const & gridShift) const
{
	assert( startBnD >= 0 );
	assert( (endBnd > startBnD) && (endBnd <= this->get_nData_gpt()) );

	LatticeStructure::RegularSymmetricGrid newGrid;
	newGrid.initialize(
			this->get_grid().interpret_fft_dim_input(newDims),
			this->get_grid().get_grid_prec(),
			gridShift,
			this->get_grid().get_symmetry(),
			this->get_grid().get_lattice());

	std::vector<double> bndData;
	std::vector<double> allData(newGrid.get_np_irred()*(endBnd-startBnD));
	for ( int ib = startBnD ; ib < endBnd; ++ib)
	{
		this->generate_interpolated_reducible_data(
				{ib},
				newGrid.view_bare_grid(),
				bndData);

		int idSymIndex = newGrid.get_symmetry().get_identity_index();
		for ( int ik = 0 ; ik < newGrid.get_np_irred(); ++ik)
			allData[ik*(endBnd-startBnD)+(ib-startBnD)] =
					bndData[newGrid.get_maps_irreducible_to_reducible()[ik][idSymIndex]];
	}

	ElectronicBands newBands;
	newBands.initialize((endBnd-startBnD), 0.0, std::move(allData), std::move(newGrid));
	return newBands;
}

void
ElectronicBands::write_tetrahedra_dos_file(
		std::string const & filename,
		std::vector<double> frequencies) const
{
	auto tetra = std::make_shared<LatticeStructure::TetrahedraGrid>();
	tetra->initialize( std::make_shared<LatticeStructure::RegularSymmetricGrid>(this->get_grid()));

	std::ofstream file( filename.c_str() );
	if ( not file.good() )
		throw std::runtime_error(std::string("Problem opening file ")+filename+" for writing the phonon DOS data");

	std::vector<double> dos;
	this->compute_DOS_tetra(tetra, frequencies, dos);

	auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
	file << "# phonon DOS in units of 1/THz with frequency in units of THz. Date is " << std::ctime(&now);

	for ( int iw = 0 ; iw < frequencies.size() ; ++iw)
		file << frequencies[iw] << '\t' << dos[iw] << '\n';
	file.close();
}

} /* namespace ElectronicStructure */
} /* namespace elephon */
