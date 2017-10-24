/*	This file NestingFunction.cpp is part of elephon.
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
 *  Created on: Oct 22, 2017
 *      Author: A. Linscheid
 */

#include "ElectronicStructure/NestingFunction.h"
#include <algorithm>
#include <assert.h>

namespace elephon
{
namespace ElectronicStructure
{


void
NestingFunction::compute_nesting_function(
		std::shared_ptr<const TetrahedraIsosurface> isoSurface)
{
	assert(isoSurface);
	isoSurface_ = isoSurface;
	nesting_.reserve( isoSurface_->get_num_iso_energies() );
	auto regularGrid = isoSurface_->get_tetra_grid()->get_grid();

	int nBnd = isoSurface_->get_nBnd();

	std::vector<int> closestGridIndices;
	std::vector<LatticeStructure::RegularSymmetricGrid::GridCube> gridCubesWithQVectors;

	double V = isoSurface_->get_tetra_grid()->get_grid()->get_lattice().get_reci_volume();
	double dV =  V / static_cast<double>(regularGrid->get_grid_dim()[0]*
										 regularGrid->get_grid_dim()[1]*
										 regularGrid->get_grid_dim()[2]);

	for (int iE = 0 ; iE < isoSurface_->get_num_iso_energies(); ++iE )
	{
		// allocate storage
		std::vector<float> nesting( nBnd*nBnd* regularGrid->get_np_red(), 0.0 );
		auto mem_layout = [&] (int ib, int ibp, int ikr)
		{
			return (ikr*nBnd+ib)*nBnd+ibp;
		};

		for (int ibnd = 0 ; ibnd < isoSurface_->get_nBnd(); ++ibnd )
		{
			std::vector<double> isoK, isoWeights;
			isoSurface_->get_irreducible_iso_vector_integration_weights_no_multiplicty(iE, ibnd, isoK, isoWeights);
			for (int ibndP = 0 ; ibndP < isoSurface_->get_nBnd(); ++ibndP )
			{
				std::vector<double> isoKP, isoWeightsP;
				isoSurface_->get_reducible_iso_vector_integration_weights(iE, ibndP, isoKP, isoWeightsP);

				for ( int ikIso = 0 ; ikIso < isoWeights.size(); ++ikIso )
				{
					// generate q = k' - k
					auto q = isoKP;
					for ( int ikIsoP = 0 ; ikIsoP < isoWeightsP.size() ; ++ikIsoP )
					{
						q[ikIsoP*3+0] -= isoK[ikIso*3+0];
						q[ikIsoP*3+1] -= isoK[ikIso*3+1];
						q[ikIsoP*3+2] -= isoK[ikIso*3+2];
					}

					// map q vectors back to the 1. BZ
					for ( auto &qi : q )
						qi -= std::floor(qi+0.5);

					// find the closes grid vectors
					regularGrid->find_closest_reducible_grid_points(q, closestGridIndices);

					for ( int iq = 0 ; iq < isoWeightsP.size() ; ++iq )
					{
						double thisWeight = isoWeights[ikIso]*isoWeightsP[iq] / dV;
						nesting[mem_layout(ibnd, ibndP, closestGridIndices[iq])] += thisWeight;
					}
				}
			}
		}

		LatticeStructure::DataRegularGrid<float> nestingGrid;
		nestingGrid.initialize_accumulation( nBnd*nBnd, 0.0, std::move(nesting), *regularGrid );
		nesting_.push_back( std::make_shared<LatticeStructure::DataRegularGrid<float>>( std::move(nestingGrid) ) );
	}
}

std::shared_ptr<LatticeStructure::DataRegularGrid<float> >
NestingFunction::get_nesting(int isoEnergyIndex) const
{
	assert( (nesting_.size() > isoEnergyIndex) and (isoEnergyIndex >= 0) );
	return nesting_[isoEnergyIndex];
}

} /* namespace ElectronicStructure */
} /* namespace elephon */
