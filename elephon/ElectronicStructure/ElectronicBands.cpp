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
#include <stdexcept>

namespace elephon
{
namespace ElectronicStructure
{

void
ElectronicBands::initialize(
		std::vector<double> const & kpoints,
		int numBands,
		std::vector<double> bandData,
		LatticeStructure::RegularGrid grid)
{
	grid_ = std::move(grid);
	assert( grid_.is_reci() );
	std::vector<int> pointIndices;
	grid_.get_list_lattice_point_indices(kpoints, pointIndices);
	nBnd_ = numBands;

	//Check for errors
	int npIrr = static_cast<int>(pointIndices.size());
	if( not ( npIrr ==  grid_.get_np_irred() ) )
		throw std::runtime_error("Initialization data must be complete in the irreducible zone.");

	if ( bandData.size() != static_cast<size_t>(numBands)*kpoints.size()/3 )
		throw std::logic_error("Called initialize inconsistent data size.");

	dataIrred_ = std::vector<double>(npIrr*numBands);
	for ( int irr = 0 ; irr < npIrr; ++irr )
		for ( int ibnd = 0 ; ibnd < nBnd_; ++ibnd )
			dataIrred_[ pointIndices[irr]*nBnd_ + ibnd ] = bandData[ irr*nBnd_ + ibnd ];
}

int
ElectronicBands::get_nBnd() const
{
	return nBnd_;
}

void
ElectronicBands::generate_reducible_grid_bands(
		std::vector<int> const & bIndices,
		std::vector<double> & bands) const
{
	std::vector<int> redIndices( grid_.get_np_red() ), irredIndices;
	for (int i = 0 ; i < grid_.get_np_red() ; ++i )
		redIndices[i] = i;

	grid_.convert_reducible_irreducible( redIndices, irredIndices );
	if ( bands.size() != bIndices.size()*redIndices.size() )
		bands = std::vector<double>(bIndices.size()*redIndices.size());
	int nBRequest = static_cast<int>(bIndices.size());
	for ( int ired = 0 ; ired < grid_.get_np_red(); ++ired )
		for ( int i = 0 ; i < nBRequest; ++i )
		{
			int ibnd = bIndices[i];
			assert( ibnd < nBnd_);
			bands[ired*nBRequest+i] = dataIrred_[ irredIndices[ired]*nBnd_ + ibnd ];
		}
}

LatticeStructure::RegularGrid const &
ElectronicBands::get_grid() const
{
	return grid_;
}

} /* namespace ElectronicStructure */
} /* namespace elephon */
