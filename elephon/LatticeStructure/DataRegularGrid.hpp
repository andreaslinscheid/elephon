/*	This file DataRegularGrid.hpp is part of elephon.
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

#include "LatticeStructure/DataRegularGrid.h"
#include "Algorithms/LocalDerivatives.h"
#include "Algorithms/FFTInterface.h"
#include <set>
#include <stdexcept>

namespace elephon
{
namespace LatticeStructure
{

template<typename T>
void
DataRegularGrid<T>::initialize(
		int numBands,
		T zeroEnergy,
		std::vector<T> bandData,
		LatticeStructure::RegularSymmetricGrid grid)
{
	grid_ = std::move(grid);
	assert( grid_.is_reci() );
	nBnd_ = numBands;

	if ( bandData.size() == grid_.get_np_irred()*nBnd_ )
	{
		// accepted as irreducible data
		dataIrred_ = std::move(bandData);
		for ( auto &e : dataIrred_ )
			e -= zeroEnergy;
	}
	else if ( bandData.size() == grid_.get_np_red()*nBnd_ )
	{
		// accepted as reducible data - convert to irreducible
		dataIrred_.resize( grid_.get_np_irred()*nBnd_ );
		for ( int irr = 0 ; irr < grid_.get_np_irred(); ++irr )
		{
			int isymId = grid_.get_symmetry().get_identity_index();
			int ireducible = grid_.get_maps_irreducible_to_reducible()[irr][isymId];
			for ( int ibnd = 0 ; ibnd < nBnd_; ++ibnd )
				dataIrred_[ irr*nBnd_ + ibnd ] = bandData[ ireducible*nBnd_ + ibnd ] - zeroEnergy;
		}
	}
	else
		throw std::runtime_error("Problem in DataRegularGrid::initialize : input data does not match grid");
}

template<typename T>
std::vector<int>
DataRegularGrid<T>::get_bands_crossing_energy_lvls(
		std::vector<double> const & energies ) const
{
	assert( energies.size() > 0 );
	std::set<int> bandset;
	for ( auto e : energies )
		for ( int ib = 0 ; ib < nBnd_; ++ib )
		{
			double refEne = double(dataIrred_[ib]) - e;
			for ( int ik = 0 ; ik < grid_.get_np_irred(); ++ik )
				if ( refEne*(dataIrred_[ik*nBnd_ + ib] - e) < 0 )
				{
					bandset.insert(ib);
					break;
				}
		}
	return std::vector<int>(bandset.begin(),bandset.end());
}

template<typename T>
std::vector<int>
DataRegularGrid<T>::get_bands_crossing_energy_window(
		std::vector<double> const & energies ) const
{
	assert( energies.size() == 2 );
	std::vector<int> bandset;
	for ( int ib = 0 ; ib < nBnd_; ++ib )
		for ( int ik = 0 ; ik < grid_.get_np_irred(); ++ik )
			if ( (dataIrred_[ik*nBnd_ + ib] >= energies[0]) and (dataIrred_[ik*nBnd_ + ib] <= energies[1]) )
			{
				bandset.push_back(ib);
				break;
			}
	return bandset;
}

template<typename T>
int
DataRegularGrid<T>::get_nBnd() const
{
	return nBnd_;
}

template<typename T>
void
DataRegularGrid<T>::generate_reducible_data(
		std::vector<int> const & bIndices,
		std::vector<T> & bands) const
{
	std::vector<int> redIndices( grid_.get_np_red() ), irredIndices;
	for (int i = 0 ; i < grid_.get_np_red() ; ++i )
		redIndices[i] = i;

	grid_.convert_reducible_irreducible( redIndices, irredIndices );
	bands.resize(bIndices.size()*redIndices.size());
	int nBRequest = static_cast<int>(bIndices.size());
	for ( int ired = 0 ; ired < grid_.get_np_red(); ++ired )
		for ( int i = 0 ; i < nBRequest; ++i )
		{
			int ibnd = bIndices[i];
			assert( ibnd < nBnd_);
			bands[ired*nBRequest+i] = dataIrred_[ irredIndices[ired]*nBnd_ + ibnd ];
		}
}

template<typename T>
void
DataRegularGrid<T>::generate_interpolated_reducible_data(
		std::vector<int> const & bIndices,
		LatticeStructure::RegularBareGrid const & interpolationGrid,
		std::vector<T> & interpolatedReducibleData) const
{
	std::vector<T> reducibleData;
	this->generate_reducible_data(bIndices, reducibleData);

	Algorithms::FFTInterface fft;
	fft.fft_interpolate(
			grid_.get_grid_dim(),
			grid_.get_grid_shift(),
			reducibleData,
			interpolationGrid.get_grid_dim(),
			interpolationGrid.get_grid_shift(),
			interpolatedReducibleData,
			bIndices.size()	);
}

template<typename T>
void
DataRegularGrid<T>::fft_interpolate(
		std::vector<int> const & newDims,
		std::vector<double> const & gridShift)
{
	assert( gridShift.size() == 3 );
	auto fftd = grid_.interpret_fft_dim_input(newDims);
	if ( (grid_.get_grid_dim() == fftd) and
			(    (std::abs(grid_.get_grid_shift()[0]-gridShift[0]) < grid_.get_grid_prec())
			 and (std::abs(grid_.get_grid_shift()[1]-gridShift[1]) < grid_.get_grid_prec())
			 and (std::abs(grid_.get_grid_shift()[2]-gridShift[2]) < grid_.get_grid_prec()) ))
		return;

	std::vector<double> oldData;
	int nB = this->get_nBnd();
	std::vector<int> bandList(nB);
	for ( int ib = 0 ; ib < nB; ++ib )
		   bandList[ib] = ib;
	this->generate_reducible_data(bandList, oldData);

	std::vector<double> newReducibleData;
	Algorithms::FFTInterface fft;
	fft.fft_interpolate(
				   grid_.get_grid_dim(),
				   grid_.get_grid_shift(),
				   oldData,
				   fftd,
				   gridShift,
				   newReducibleData,
				   nB);

	LatticeStructure::RegularSymmetricGrid newGrid;
	newGrid.initialize(
				   fftd,
				   grid_.get_grid_prec(),
				   gridShift,
				   grid_.get_symmetry(),
				   grid_.get_lattice());

	this->initialize(
				   nB,
				   0.0,
				   newReducibleData,
				   newGrid);
}

template<typename T>
template<typename TD>
void
DataRegularGrid<T>::compute_derivatives_sqr_polynom(
		std::vector<int> const & bandIndices,
		std::vector<int> const & reducibleKPTIndices,
		std::vector<TD> * gradientFieldPtr,
		std::vector<TD> * hessianFieldPtr ) const
{
	// this lambda take the job of mapping a regular reducible k grid index and a "local" band index
	// into the irreducible zone and the full band-context band index.
	auto translate_bands = [&] (int ikr, int ib) {
		return this->read( grid_.get_maps_red_to_irreducible()[ikr], bandIndices[ib]);
	};

	Algorithms::localDerivatives::compute_derivatives_sqr_polynom<TD>(
			bandIndices.size(),
			reducibleKPTIndices,
			gradientFieldPtr,
			hessianFieldPtr,
			grid_.view_bare_grid(),
			translate_bands );
}

template<typename T>
LatticeStructure::RegularSymmetricGrid const &
DataRegularGrid<T>::get_grid() const
{
	return grid_;
}

template<typename T>
std::pair<T, T>
DataRegularGrid<T>::get_min_max() const
{
	if ( dataIrred_.empty() )
		return std::make_pair(0.0, 0.0);
	auto mm = std::minmax_element(dataIrred_.begin(), dataIrred_.end());

	return std::make_pair(*mm.first, *mm.second);
}

template<typename T>
T
DataRegularGrid<T>::read(int i, int ib) const
{
	assert( (ib >= 0) && (ib < nBnd_));
	assert( (i*nBnd_+ib >= 0) && (i*nBnd_+ib < dataIrred_.size()));
	return dataIrred_[i*nBnd_+ib];
}

template<typename T>
T &
DataRegularGrid<T>::write(int i, int ib)
{
	assert( (ib >= 0) && (ib < nBnd_));
	assert( (i*nBnd_+ib >= 0) && (i*nBnd_+ib < dataIrred_.size()));
	return dataIrred_[i*nBnd_+ib];
}

} /* namespace LatticeStructure */
} /* namespace elephon */
