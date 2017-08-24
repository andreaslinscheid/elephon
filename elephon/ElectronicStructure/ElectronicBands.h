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

#include <vector>
#include "LatticeStructure/RegularSymmetricGrid.h"

namespace elephon
{
namespace ElectronicStructure
{

class ElectronicBands
{
public:

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
	 * Set the data in the the band structure.
	 *
	 * The data is assumed to match the shape of either the reducible or the irreducible grid.
	 * Which one is picked depends on the size of \p bandData ; if it matches either one it will be chosen.
	 *
	 * @param numBands	Number of bands
	 * @param fermiEnergy	energy reference of the \p bandData
	 * @param bandData	N, which is the number of k points, times number of bands data values
	 * 					with 'bands' as the fast running dimension such as [(b=0,k=0),(b=1,k=0),...(b=num bands,k=N-1)]
	 * @param grid		K point grid
	 */
	void initialize(
			int numBands,
			double fermiEnergy,
			std::vector<double> bandData,
			LatticeStructure::RegularSymmetricGrid grid);

	std::vector<int> get_bands_crossing_energy_lvls(
			std::vector<double> const & energies ) const;

	int get_nBnd() const;

	int get_nspin() const;

	void generate_reducible_grid_bands(
			std::vector<int> const & bIndices,
			std::vector<double> & bands) const;

	LatticeStructure::RegularSymmetricGrid const & get_grid() const;

	/**
	 * Retrieve a band data value at a point in the grid.
	 *
	 * @param ikIrred	The irreducible index of the k grid point.
	 * @param ib		The band index [0,nBnd[
	 * @param ispin		The spin index
	 * @return	Band energy value in eV relative to the Fermi level at E=0.
	 */
	double operator() (int ikIrred, int ib, int ispin = 0) const;

	void fft_interpolate(std::vector<int> const & newDims, std::vector<double> const & gridShift);
private:

	int nBnd_;

	std::vector<double> dataIrred_;

	LatticeStructure::RegularSymmetricGrid grid_;
};

} /* namespace ElectronicStructure */
} /* namespace elephon */

#endif /* ELEPHON_ELECTRONICSTRUCTURE_ELECTRONICBANDS_H_ */
