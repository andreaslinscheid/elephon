/*	This file ElectronicBands.h is part of elephon.
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

#ifndef ELEPHON_ELECTRONICSTRUCTURE_ELECTRONICBANDS_H_
#define ELEPHON_ELECTRONICSTRUCTURE_ELECTRONICBANDS_H_

#include "LatticeStructure/DataRegularGrid.h"
#include "LatticeStructure/RegularSymmetricGrid.h"
#include <vector>

namespace elephon
{
namespace ElectronicStructure
{

class ElectronicBands : public LatticeStructure::DataRegularGrid<double>
{
public:

	using LatticeStructure::DataRegularGrid<double>::initialize;

	/**
	 * Set the data in the the band structure.
	 *
	 * The data is located by the k point coordinates.
	 *
	 * @param kpoints	K point in direct coordinates as [k0_x,k0_y,k0_z, k1_x,...,kN-1_z]
	 * @param numBands	Number of bands
	 * @param fermiEnergy	energy reference of the \p bandData
	 * @param bandData	N, which is the number of k points, times number of bands data values
	 * 					with 'bands' as the fast running dimension such as [(b=0,k=0),(b=1,k=0),...(b=num bands,k=N-1)]
	 * @param grid		K point grid
	 */
	void initialize(
			std::vector<double> const & kpoints,
			int numBands,
			double fermiEnergy,
			std::vector<double> bandData,
			LatticeStructure::RegularSymmetricGrid grid);

	/**
	 * Retrieve the number of spin channels.
	 *
	 * @return integer with number of spin channels.
	 */
	int get_nspin() const;

	/**
	 * Retrieve the number of bands.
	 *
	 * @return integer with number of bands.
	 */
	int get_nBnd() const;

	/**
	 * Retrieve a band data value at a point in the grid.
	 *
	 * @param ikIrred	The irreducible index of the k grid point.
	 * @param ib		The band index [0,nBnd[
	 * @param ispin		The spin index
	 * @return	Band energy value in eV relative to the Fermi level at E=0.
	 */
	double operator() (int ikIrred, int ib, int ispin = 0) const;

	/**
	 * Generate a new object with FFT interpolated data of selected bands.
	 *
	 * @param startBnD	relative new first band. Must be >= 0.
	 * @param endBnd	one past the relative new last band. Must be > \p startBnD and <= this->get_nBnd() .
	 * @param newDims	1) A list of 3 positive numbers representing the new grid dimension. If a dimension is '0'
	 * 						the number will be replaced by the internal grid number (i.e. remains unchanged).
	 * 					2) A single positive number in which case the internal grid will be scaled in each direction.
	 * @param gridShift	A list of 3 floating values in the range [0,1[ representing a shift within a single cube of the
	 * 					grid. Thus, to know the shift in absolute numbers, multiply by the inverse grid dimension e.g.
	 * 					for a grid 5 5 5 with shift 0.5 0.0 0.0 the absolute shift will be (1.0/5.0)*0.5 in x direction.
	 * @return A new band object with interpolated data.
	 */
	ElectronicBands fft_interpolate_part(
			int startBnD, int endBnd,
			std::vector<int> const & newDims,
			std::vector<double> const & gridShift) const;

	/**
	 * Compute and write the electron DOS to file.
	 *
	 * @param filename		The name of the file. Will be overwritten if exists and
	 * @param energySamples	Energies where the DOS will be evaluated.
	 */
	void write_tetrahedra_dos_file(
			std::string const & filename,
			std::vector<double> energySamples) const;
};

} /* namespace ElectronicStructure */
} /* namespace elephon */

#include "ElectronicStructure/ElectronicBands.hpp"
#endif /* ELEPHON_ELECTRONICSTRUCTURE_ELECTRONICBANDS_H_ */
