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
#include <stdexcept>
#include <set>

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
	grid_ = std::move(grid);
	assert( grid_.is_reci() );
	std::vector<int> pointIndices;
	grid_.get_list_lattice_point_indices(kpoints, pointIndices);
	nBnd_ = numBands;

	//Check for errors
	int npIrr = pointIndices.size();
	if( not ( npIrr ==  grid_.get_np_irred() ) )
		throw std::runtime_error("Initialization data must be complete in the irreducible zone.");

	if ( bandData.size() != static_cast<size_t>(numBands)*kpoints.size()/3 )
		throw std::logic_error("Called initialize inconsistent data size.");

	dataIrred_ = std::vector<double>(npIrr*numBands);
	for ( int irr = 0 ; irr < npIrr; ++irr )
		for ( int ibnd = 0 ; ibnd < nBnd_; ++ibnd )
			dataIrred_[ pointIndices[irr]*nBnd_ + ibnd ] = bandData[ irr*nBnd_ + ibnd ] - fermiEnergy;
}

void
ElectronicBands::initialize(
		int numBands,
		double fermiEnergy,
		std::vector<double> bandData,
		LatticeStructure::RegularSymmetricGrid grid)
{
	grid_ = std::move(grid);
	assert( grid_.is_reci() );
	nBnd_ = numBands;

	if ( bandData.size() == grid_.get_np_irred()*nBnd_ )
	{
		dataIrred_ = std::move(bandData);
		for ( auto &e : dataIrred_ )
			e -= fermiEnergy;
	}
	else if ( bandData.size() == grid_.get_np_red()*nBnd_ )
	{
		dataIrred_.resize( grid_.get_np_irred()*nBnd_ );
		for ( int ired = 0 ; ired < grid_.get_np_red(); ++ired )
		{
			int irr = grid_.get_maps_red_to_irreducible()[ired];
			for ( int ibnd = 0 ; ibnd < nBnd_; ++ibnd )
			{
				dataIrred_[ irr*nBnd_ + ibnd ] = bandData[ ired*nBnd_ + ibnd ] - fermiEnergy;
			}
		}
	}
	else
		throw std::runtime_error("Problem in ElectronicBands::initialize : input data does not match grid");
}

int
ElectronicBands::get_nBnd() const
{
	return nBnd_;
}

int
ElectronicBands::get_nspin() const
{
	return 1;
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

std::vector<int>
ElectronicBands::get_bands_crossing_energy_lvls(
		std::vector<double> const & energies ) const
{
	assert( energies.size() > 0 );
	std::set<int> bandset;
	for ( auto e : energies )
		for ( int ib = 0 ; ib < nBnd_; ++ib )
		{
			double refEne = dataIrred_[ib] - e;
			for ( int ik = 0 ; ik < grid_.get_np_irred(); ++ik )
				if ( refEne*(dataIrred_[ik*nBnd_ + ib] - e) < 0 )
				{
					bandset.insert(ib);
					break;
				}
		}
	return std::vector<int>(bandset.begin(),bandset.end());
}

LatticeStructure::RegularSymmetricGrid const &
ElectronicBands::get_grid() const
{
	return grid_;
}

double
ElectronicBands::operator() (int ikIrred, int ib, int ispin) const
{
	assert(ikIrred*nBnd_ + ib < dataIrred_.size());
	return dataIrred_[ikIrred*nBnd_ + ib];
}

void
ElectronicBands::fft_interpolate(
		std::vector<int> const & newDims,
		std::vector<double> const & gridShift)
{
       auto fftd = newDims;
       if ( fftd.size() == 1 )
       {
               int scale = fftd.at(0);
               fftd = grid_.get_grid_dim();
               if ( scale != 0 )
                       for ( auto &d : fftd )
                               d *= scale;
       }
       else
       {
               if ( newDims.size() != 3 )
                       throw std::runtime_error("Incorrect grid dimension for bands interpolate.");
               for ( int id = 0 ; id < 3; ++id)
                       fftd[id] = fftd[id] == 0 ? grid_.get_grid_dim()[id] : fftd[id];
       }

       std::vector<double> oldData;
       int nB = this->get_nBnd();
       std::vector<int> bandList(nB);
       for ( int ib = 0 ; ib < nB; ++ib )
               bandList[ib] = ib;
       this->generate_reducible_grid_bands(bandList, oldData);

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


} /* namespace ElectronicStructure */
} /* namespace elephon */
