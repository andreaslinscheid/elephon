/*	This file GridRotationMap.cpp is part of elephon.
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
 *  Created on: Sep 26, 2017
 *      Author: A. Linscheid
 */

#include "Algorithms/GridRotationMap.h"

namespace elephon
{
namespace Algorithms
{

void compute_grid_rotation_map(
		std::vector<double> const & shift,
		LatticeStructure::RegularBareGrid const & grid,
		LatticeStructure::Symmetry const & symmetry,
		std::vector< std::vector<int> > & rotMap)
{
	int nR = grid.get_num_points();
	int iS = symmetry.get_num_symmetries();
	rotMap = std::vector< std::vector<int> >(iS, std::vector<int>(nR) );
	std::vector<double> shiftedGrid(nR*3);
	for ( int i = 0 ; i < nR ; ++i)
	{
		auto r = grid.get_vector_direct(i);
		for ( int xi = 0 ; xi < 3; ++xi)
			shiftedGrid[i*3+xi] = r[xi] - shift[xi];
	}

	for ( int isym = 0 ; isym < iS; ++isym)
	{
		std::vector<double> shiftedGridCpy = shiftedGrid;
		//rotate all grid points
		symmetry.apply(isym, shiftedGridCpy.begin(), shiftedGridCpy.end(), true);
		//since the grid is regular, we must obtain the (x,y,z) indices by scaling with the grid dimension
		//in each direction
		std::vector<int> xyz(3);
		for ( int ir = 0 ; ir < nR; ++ir)
		{
			for ( int j = 0 ; j < 3; ++j)
			{
				shiftedGridCpy[ir*3+j] -= std::floor(shiftedGridCpy[ir*3+j]);
				shiftedGridCpy[ir*3+j] *= grid.get_grid_dim()[j];
				xyz[j] = std::floor(shiftedGridCpy[ir*3+j]+0.5);
				if ( std::abs(shiftedGridCpy[ir*3+j]-xyz[j]) > 0.01 )
					throw std::logic_error("The symmetry operation does not map the grid to itself which can't be.");
				xyz[j] = xyz[j] >= grid.get_grid_dim()[j] ? xyz[j] - grid.get_grid_dim()[j] : xyz[j];
			}
			int cnsq = grid.get_xyz_to_reducible(xyz);
			assert( (cnsq >= 0) and (cnsq < nR));
			rotMap[isym][cnsq] = ir;
		}
	}
}

} /* namespace Algorithms */
} /* namespace elephon */