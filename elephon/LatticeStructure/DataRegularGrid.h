/*	This file DataRegularGrid.h is part of elephon.
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
 *  Created on: Oct 2, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_LATTICESTRUCTURE_DATAREGULARGRID_H_
#define ELEPHON_LATTICESTRUCTURE_DATAREGULARGRID_H_

#include "LatticeStructure/RegularSymmetricGrid.h"
#include <vector>

namespace elephon
{
namespace LatticeStructure
{

template<typename T>
class DataRegularGrid
{
public:
	/**
	 * Set the data in the the data structure.
	 *
	 * The data is assumed to match the shape of either the reducible or the irreducible grid.
	 * Which one is picked depends on the size of \p bandData ; if it matches either one it will be chosen.
	 *
	 * @param numBands	Number of bands
	 * @param zeroEnergy	energy reference of the \p bandData
	 * @param bandData	N, which is the number of k points, times number of bands data values
	 * 					with 'bands' as the fast running dimension such as [(b=0,k=0),(b=1,k=0),...(b=num bands,k=N-1)]
	 * @param grid		K point grid
	 */
	void initialize(
			int numBands,
			T zeroEnergy,
			std::vector<T> bandData,
			LatticeStructure::RegularSymmetricGrid grid);

	std::vector<int> get_bands_crossing_energy_lvls(
			std::vector<double> const & energies ) const;

	std::vector<int> get_bands_crossing_energy_window(
			std::vector<double> const & energies ) const;

	int get_nBnd() const;

	void generate_reducible_data(
			std::vector<int> const & bIndices,
			std::vector<T> & bands) const;

	void generate_interpolated_reducible_data(
			std::vector<int> const & bIndices,
			LatticeStructure::RegularBareGrid const & interpolationGrid,
			std::vector<T> & interpolatedReducibleData) const;


	LatticeStructure::RegularSymmetricGrid const & get_grid() const;

	/**
	 * Replace the content of this object with FFT interpolated data.
	 *
	 * If the grids are determined to be the same, no action is taken.
	 *
	 * @param newDims	1) A list of 3 positive numbers representing the new grid dimension. If a dimension is '0'
	 * 						the number will be replaced by the internal grid number (i.e. remains unchanged).
	 * 					2) A single positive number in which case the internal grid will be scaled in each direction.
	 * @param gridShift	A list of 3 floating values in the range [0,1[ representing a shift within a single cube of the
	 * 					grid. Thus, to know the shift in absolute numbers, multiply by the inverse grid dimension e.g.
	 * 					for a grid 5 5 5 with shift 0.5 0.0 0.0 the absolute shift will be (1.0/5.0)*0.5 in x direction.
	 */
	void fft_interpolate(
			std::vector<int> const & newDims,
			std::vector<double> const & gridShift);

	/**
	 * Compute the minimal and maximal value of the band data.
	 *
	 * @return	pair of first minimal and second maximal value in eV
	 */
	std::pair<T, T> get_min_max() const;

	template<typename TD>
	void compute_derivatives_sqr_polynom(
			std::vector<int> const & bandIndices,
			std::vector<int> const & reducibleKPTIndices,
			std::vector<TD> * gradientFieldPtr,
			std::vector<TD> * hessianFieldPtr ) const;

	T read(int i, int ib) const;

	T & write(int i, int ib);
private:

	int nBnd_ = 0;

	std::vector<T> dataIrred_;

	LatticeStructure::RegularSymmetricGrid grid_;
};

} /* namespace LatticeStructure */
} /* namespace elephon */

#include "LatticeStructure/DataRegularGrid.hpp"
#endif /* ELEPHON_LATTICESTRUCTURE_DATAREGULARGRID_H_ */
