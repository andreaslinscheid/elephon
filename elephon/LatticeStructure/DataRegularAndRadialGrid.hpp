/*	This file DataRegularAndRadialGrid.hpp is part of elephon.
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

#include "LatticeStructure/DataRegularAndRadialGrid.h"

namespace elephon
{
namespace LatticeStructure
{

template<typename T>
void
DataRegularAndRadialGrid<T>::initialize(
		LatticeStructure::RegularBareGrid regularGrid,
		Auxillary::alignedvector::aligned_vector<T> regularGridData,
		std::vector<AtomicSite::AtomSiteData> radialGridData)
{
	regularGrid_ = std::move(regularGrid);
	regularGridData_ = std::move(regularGridData);
	radialGridData_ = std::move(radialGridData);
}

template<typename T>
int
DataRegularAndRadialGrid<T>::get_max_num_radial_elements() const
{
	int nRadMax = 0;
	for (auto const & dv : radialGridData_)
		nRadMax = std::max(dv.get_data().get_radial_grid().get_num_R(),nRadMax);
	return nRadMax;
}

template<typename T>
int
DataRegularAndRadialGrid<T>::get_max_num_angular_moment_channels() const
{
	auto lM = this->get_max_angular_moment();
	return Auxillary::memlayout::angular_momentum_layout(lM+1,-lM-1);;
}

template<typename T>
int
DataRegularAndRadialGrid<T>::get_max_angular_moment() const
{
	int lmax = 0;
	for (auto const & dv : radialGridData_)
		lmax = std::max(lmax, dv.get_data().get_l_max());
	return lmax;
}

template<typename T>
void
DataRegularAndRadialGrid<T>::transform(symmetry::SymmetryOperation const & sop)
{
	sop.transform_scalar_field_regular_grid(regularGrid_, regularGridData_);
	// apply to each atom data, then reshuffle the atoms
	for (auto & ad : radialGridData_)
		ad.edit_data().transform(sop);

	// todo this logic below is already implemented in the atom symmetry mappings
	// unfortunately its not easy to get this information here without changing the logic of how
	// a symmetry operation class works.
	std::map<LatticeStructure::Atom,int> lookup;
	for (int ia = 0 ; ia < radialGridData_.size(); ++ia)
		lookup.insert(std::make_pair(radialGridData_[ia].get_atom(), ia));

	std::vector<int> symMap(radialGridData_.size()); // tells where atom index i was before applying the symmetry operation
	for (int ia = 0 ; ia < radialGridData_.size(); ++ia)
	{
		auto a = radialGridData_[ia].get_atom();
		a.transform(sop);
		auto ret = lookup.find(a);
		if (ret == lookup.end())
			throw std::logic_error(" awkwardly doubled atom symmetry mapping failed");
		symMap[ia] = ret->second;
	}

	// since the atoms are mapped into each other and we keep the order, only the data
	// has to be moved around
	decltype(radialGridData_) radialGridDataBuffer(radialGridData_.size());
	for (int ia = 0 ; ia < radialGridData_.size(); ++ia)
	{
		const int iaMapped = symMap[ia];
		radialGridDataBuffer[ia].initialize(radialGridData_[ia].get_atom(), std::move(radialGridData_[iaMapped].edit_data()));
	}

	std::swap(radialGridData_,radialGridDataBuffer);
}

template<typename T>
typename Auxillary::alignedvector::aligned_vector<T>::const_iterator
DataRegularAndRadialGrid<T>::begin_regular_data() const
{
	return regularGridData_.begin();
}

template<typename T>
typename Auxillary::alignedvector::aligned_vector<T>::const_iterator
DataRegularAndRadialGrid<T>::end_regular_data() const
{
	return regularGridData_.end();
}

template<typename T>
Auxillary::alignedvector::ZV::const_iterator
DataRegularAndRadialGrid<T>::begin_radial_data(int atomIndex) const
{
	assert((atomIndex>=0)&&(atomIndex<radialGridData_.size()));
	return radialGridData_[atomIndex].get_data().begin();
}

template<typename T>
Auxillary::alignedvector::ZV::const_iterator
DataRegularAndRadialGrid<T>::end_radial_data(int atomIndex) const
{
	assert((atomIndex>=0)&&(atomIndex<radialGridData_.size()));
	return radialGridData_[atomIndex].get_data().end();
}

template<typename T>
int
DataRegularAndRadialGrid<T>::get_num_atom_data_sets() const
{
	return radialGridData_.size();
}

template<typename T>
AtomicSite::AtomSiteData const &
DataRegularAndRadialGrid<T>::view_radial_data_set(int index) const
{
	assert((index >=0)&&(index < this->get_num_atom_data_sets()));
	return radialGridData_[index];
}

} /* namespace LatticeStructure */
} /* namespace elephon */
