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
#include "LatticeStructure/TetrahedraGrid.h"
#include "Auxillary/AlignedVector.h"
#include <vector>
#include <memory>

namespace elephon
{
namespace LatticeStructure
{

template<typename T>
class DataRegularGrid
{
public:

	/**
	 * Define a generic vector type for memory aligned data.
	 */
	typedef Auxillary::alignedvector::aligned_vector<T,Auxillary::alignedvector::architecture_align> VT;

	typedef Auxillary::alignedvector::aligned_vector<std::complex<T>,Auxillary::alignedvector::architecture_align> VCT;

	template<class F>
	void initialize(
			T referenceEnergy,
			F const & functor,
			LatticeStructure::RegularSymmetricGrid grid);

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
			VT bandData,
			LatticeStructure::RegularSymmetricGrid grid);

	/**
	 * Set the data in the the data structure by adding data into the irreducible zone.
	 *
	 * The data is assumed to match the reducible irreducible grid.
	 * As opposed to ::initialize it will assume the data is distributed in the reducible zone and needs to be
	 * summed. Note that \p zeroEnergy is subtracted _after_ the accumulation.
	 *
	 * @param numBands	Number of bands
	 * @param zeroEnergy	energy reference of the \p bandData
	 * @param bandData	N, which is the number of k points, times number of bands data values
	 * 					with 'bands' as the fast running dimension such as [(b=0,k=0),(b=1,k=0),...(b=num bands,k=N-1)]
	 * @param grid		K point grid
	 */
	void initialize_accumulation(
			int numBands,
			T zeroEnergy,
			VT bandData,
			LatticeStructure::RegularSymmetricGrid grid);

	/**
	 * return a list of bands that cross either one of submitted energy levels.
	 *
	 * @param energies	List of energy levels. Must not be empty.
	 * @return			vector of bands that cross one of these energy levels.
	 */
	std::vector<int> get_bands_crossing_energy_lvls(
			std::vector<double> const & energies ) const;

	/**
	 * return a list of bands that cross the range of the energy window.
	 *
	 * @param energies	A list of two values determining a non-empty range in energy.
	 * @return			vector of bands that cross into the energy window.
	 */
	std::vector<int> get_bands_crossing_energy_window(
			std::vector<double> const & energies ) const;

	/**
	 * return the number of data points per grind point.
	 *
	 * @return the number of points per grid point.
	 */
	int get_nData_gpt() const;

	/**
	 * compute data for bIndices in the reducible grid.
	 *
	 * @param bIndices	The band indices for which the data will be mapped onto the reducible grid.
	 * @param bands		Data for each band in the list for the reducible grid with band-major order in some form of vector VR.
	 */
	template<typename VR>
	void generate_reducible_data(
			std::vector<int> const & bIndices,
			VR & bands) const;

	/**
	 * compute data for bIndices and interpolate in to the reducible grid.
	 *
	 * @param bIndices					The band indices for which the data will be mapped onto the reducible grid.
	 * @param interpolationGrid			The reducible grid onto which the data will be interpolated.
	 * @param interpolatedReducibleData	Data for each band in the list for the reducible grid with band-major order.
	 */
	void generate_interpolated_reducible_data(
			std::vector<int> const & bIndices,
			LatticeStructure::RegularBareGrid const & interpolationGrid,
			VT & interpolatedReducibleData) const;

	/**
	 * Constant access to the internal copy of the regular grid.
	 *
	 * @return	constant reference to the regular grid.
	 */
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

	/**
	 * Compute local derivatives of the data in the grid.
	 *
	 * @param bandIndices			Computed at each grid index for these band indices.
	 * @param reducibleKPTIndices	Computed for each band index at these grid indices.
	 * @param gradientFieldPtr		If the pointer is not Null, the local gradient will be placed in the vector pointed to
	 * 								by this pointer. Reallocated as needed. Data will be in the layout dx,dy,dz for, first
	 * 								each band index and then for each grid index.
	 * @param hessianFieldPtr		If the pointer is not Null, the local hessian will be placed in the vector pointed to
	 * 								by this pointer. Reallocated as needed. Data will be in the layout
	 * 								dx*dx, dx*dy, dx*dz, dy*dy, dy*dz, dz*dz for, first
	 * 								each band index and then for each grid index. We imply the symmetry between dx*dy and dy*dx
	 */
	template<typename TD>
	void compute_derivatives_sqr_polynom(
			std::vector<int> const & bandIndices,
			std::vector<int> const & reducibleKPTIndices,
			std::vector<TD> * gradientFieldPtr,
			std::vector<TD> * hessianFieldPtr ) const;

	/**
	 * Access the internal data
	 *
	 * @param i		irreducible grid point index
	 * @param ib	band index
	 * @return		the value at this grid point / band index pair
	 */
	T read(int i, int ib) const;

	/**
	 * Access the internal data
	 *
	 * @param i		irreducible grid point index
	 * @param ib	band index
	 * @return		mutable reference to the value at this grid point / band index pair
	 */
	T & write(int i, int ib);

	std::vector<double> setup_frequency_grid(std::vector<double> range, int numpts) const;

	/**
	 * Mechanism to convert the input option 'range' into an actual energy range.
	 *
	 * Since one possible "default" is to use all the range that is there, this
	 * input parameter interpretation must be closely tied to the data, thus it
	 * is implemented here
	 *
	 * @param range	vector of no more than 2 elements according to the specifications in
	 * 				IOMethods::InputOptions for the input parameter 'ewinbnd' and 'phrange'
	 * @return	pair with the first being the minium and the second the maxium of the range.
	 */
	std::pair<T,T> interpret_range(std::vector<double> range) const;

	void compute_DOS(
			std::vector<double> const & energies,
			std::vector<T> & dos) const;

	void compute_DOS_tetra(
			std::shared_ptr<const TetrahedraGrid> tetraGrid,
			std::vector<double> const & energies,
			std::vector<T> & dos) const;

	template<class F>
	void compute_DOS_wan(
			F const & functor,
			std::vector<double> const & energies,
			std::vector<T> & dos) const;

	void interpolate_bands_along_path(
			std::vector<double> const & nonGridPoints,
			std::vector<T> energyRange,
			VT & bands,
			int &numBands,
			std::shared_ptr<const LatticeStructure::TetrahedraGrid> interpolMesh = nullptr) const;
private:

	int nDGP_ = 0;

	VT dataIrred_;

	LatticeStructure::RegularSymmetricGrid grid_;

	mutable std::shared_ptr<LatticeStructure::TetrahedraGrid> tetraGrid_;

	template<class F>
	void compute_DOS_general(
			F const & functor,
			std::vector<double> const & energies,
			std::vector<T> & dos) const;
};

} /* namespace LatticeStructure */
} /* namespace elephon */

#include "LatticeStructure/DataRegularGrid.hpp"
#endif /* ELEPHON_LATTICESTRUCTURE_DATAREGULARGRID_H_ */
