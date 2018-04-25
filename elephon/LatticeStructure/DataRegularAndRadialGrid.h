/*	This file DataRegularAndRadialGrid.h is part of elephon.
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
 *  Created on: Feb 12, 2018
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_LATTICESTRUCTURE_DATAREGULARANDRADIALGRID_H_
#define ELEPHON_LATTICESTRUCTURE_DATAREGULARANDRADIALGRID_H_

#include "Auxillary/AlignedVector.h"
#include "AtomicSite/AtomSiteData.h"
#include "symmetry/SymmetryOperation.h"
#include "LatticeStructure/RegularBareGrid.h"

namespace elephon
{

namespace LatticeStructure
{

/**
 * A container class for both regular grid and Radial grid data.
 */
template<typename T>
class DataRegularAndRadialGrid
{
public:

	/**
	 * Set this object by providing the data explicitly.
	 *
	 * @param[in] regularGrid		The grid on which the regular grid lives.
	 * @param[in] regularGridData	Copy of the regular grid data.
	 * @param[in] radialGridData	Copy of the radial grid data.
	 */
	void initialize(
			LatticeStructure::RegularBareGrid regularGrid,
			Auxillary::alignedvector::aligned_vector<T> regularGridData,
			std::vector<AtomicSite::AtomSiteData> radialGridData);

	/**
	 * Get max number of radial data elements stored in the object.
	 * @return	the number of radial data elements
	 */
	int get_max_num_radial_elements() const;

	/**
	 * Get max number of angular momentum channels of the data stored in the object.
	 * @return	the number of angular momentum channels
	 */
	int get_max_num_angular_moment_channels() const;

	/**
	 * Get max angular moment
	 * @return	The maximal angular moment quantum number.
	 */
	int get_max_angular_moment() const;

	/**
	 * Re-shuffel the internal data according to the symmetry operation sop.
	 *
	 * @param[in] sop	Representation of the symmetry operation that will be applied.
	 */
	void transform(symmetry::SymmetryOperation const & sop);

	typename Auxillary::alignedvector::aligned_vector<T>::const_iterator begin_regular_data() const;

	typename Auxillary::alignedvector::aligned_vector<T>::const_iterator end_regular_data() const;

	Auxillary::alignedvector::ZV::const_iterator begin_radial_data(int atomIndex) const;

	Auxillary::alignedvector::ZV::const_iterator end_radial_data(int atomIndex) const;

	/**
	 * Get number of angular data sets.
	 * @return	number of atomic site data sets.
	 */
	int get_num_atom_data_sets() const;

	/**
	 * Inspect a radial data set stored in this object.
	 * @param[in] index 	index in the list of atoms, i.e. data sets.
	 * @return	A const reference to the AtomicSite::AtomSiteData
	 */
	AtomicSite::AtomSiteData const & view_radial_data_set(int index) const;
private:

	int numElementsRadialTotal_ = 0;

	LatticeStructure::RegularBareGrid regularGrid_;

	Auxillary::alignedvector::aligned_vector<T> regularGridData_;

	std::vector<AtomicSite::AtomSiteData> radialGridData_;
};

} /* namespace LatticeStructure */
} /* namespace elephon */

#include "LatticeStructure/DataRegularAndRadialGrid.hpp"
#endif /* ELEPHON_LATTICESTRUCTURE_DATAREGULARANDRADIALGRID_H_ */
