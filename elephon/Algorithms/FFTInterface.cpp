/*	This file FFTInterface.cpp is part of elephon.
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
 *  Created on: Jul 7, 2017
 *      Author: A. Linscheid
 */

#include "Algorithms/FFTInterface.h"
#include <assert.h>

namespace elephon
{
namespace Algorithms
{

FFTInterface::FFTInterface()
{

}

void
FFTInterface::plan_fft(
		std::vector<int> const & gridDims,
		int exponentSign,
		int nDataPerGridPt,
		int hintHowOften,
		fftw_plan_s * & toAlloc,
		bool dataLayoutRowMajor)
{
	int gridnum = 1;
	for ( auto nd : gridDims )
		gridnum *= nd;

	//Please see the fftw3 documentation for the meaning of the variables below.
	//I am using the same nomenclature.
	//TODO This can be done slightly better using the real to complex FFT
	int dim = gridDims.size();
	int * n = new int [ dim ];
	for ( int i = 0 ; i < dim ; ++i )
		n[i] = dataLayoutRowMajor ? gridDims[i] : gridDims[dim - 1 - i] ;
	int dblock = gridnum;
	int * inembed =NULL;
	int * onembed = NULL;
	int howmany = nDataPerGridPt;
	int idist = dblock;
	int odist = dblock;
	int istride = 1;
	int ostride = 1;

	auto plan_how_many = hintHowOften > 2 ? (hintHowOften > 10 ? FFTW_PATIENT : FFTW_MEASURE) : FFTW_ESTIMATE;

	toAlloc = fftw_plan_many_dft(
				dim,n,howmany,
				FFTBuffer_, inembed,istride, idist,
				FFTBuffer_, onembed,ostride, odist,
				exponentSign,plan_how_many);

	delete [] n;
	delete [] inembed;
	delete [] onembed;
}

FFTInterface::~FFTInterface()
{
	fftw_free( FFTBuffer_ );
	fftw_destroy_plan( fftw3PlanFowdRtoK_ );
	fftw_destroy_plan( fftw3PlanBkwdKtoR_ );
}

void
FFTInterface::inplace_to_freq(int &x, int &y, int &z, int dx, int dy, int dz)
{
	assert((x >= 0 && x < dx) && (y >= 0 && y < dy) && (z >= 0 && z < dz));
	x = x <= dx/2 ? x : x - dx;
	y = y <= dy/2 ? y : y - dy;
	z = z <= dz/2 ? z : z - dz;
}

void
FFTInterface::inplace_to_freq(std::vector<int> & g, std::vector<int> const & d)
{
	assert(g.size() == d.size());
	for ( int i = 0 ; i < g.size(); ++i )
	{
		assert((g[i] >= 0 ) && (g[i] < d[i]));
		g[i] = g[i] <= d[i]/2 ? g[i] : g[i] - d[i];
	}
}

void
FFTInterface::freq_to_inplace(int &x, int &y, int &z, int dx, int dy, int dz)
{
	assert((x > -dx/2-dx%2 && x <= dx/2) && (y > -dy/2-dy%2 && y <= dy/2) && (z > -dz/2-dz%2 && z <= dz/2));
	x = x < 0 ? x + dx : x;
	y = y < 0 ? y + dy : y;
	z = z < 0 ? z + dz : z;
}

void
FFTInterface::freq_to_inplace(std::vector<int> & g, std::vector<int> const & d)
{
	assert(g.size() == d.size());
	for ( int i = 0 ; i < g.size(); ++i )
	{
		assert((g[i] > -d[i]/2-d[i]%2) && (g[i] <= d[i]/2));
		g[i] = g[i] < 0 ? g[i] + d[i]: g[i];
	}
}

void
FFTInterface::cnsq_to_xyz(
		int cnsq,
		std::vector<int> & xyz,
		std::vector<int> const& grid,
		bool dataLayoutRowMajor)
{
	int d = grid.size();
	assert(xyz.size()==d);
	std::vector<int> dimProd( d , 1 );
	std::fill(xyz.begin(), xyz.end(), cnsq );
	if (not dataLayoutRowMajor)
	{
		//reverse-engineer expressions like
		//	index = x + dx*( y + dy*( z ) )
		//	      = x + dx*y + dx*dy*z
		for (int i = 1 ; i < d ; ++i)
			dimProd[i] = grid[ i-1 ]*dimProd[i-1];

		for ( int i = d-1 ; i >= 1; --i )
		{
			xyz[i] = xyz[i]/dimProd[i];
			for ( int j = i-1 ; j >=0; --j )
				xyz[j] -= xyz[i]*dimProd[i];
		}
	}
	else
	{
		//reverse-engineer expressions like
		//	index = (x*dy + y)*dz + z
		//		  = x*dy*dz + y*dz + z
		for (int i = d-2 ; i >= 0 ; --i)
			dimProd[i] = grid[ i+1 ]*dimProd[i+1];

		for ( int i = 0 ;  i < d-1 ; ++i )
		{
			xyz[i] = xyz[i]/dimProd[i];
			for ( int j = i+1 ;  j < d ; ++j )
				xyz[j] -= xyz[i]*dimProd[i];
		}
	}
}

void
FFTInterface::xyz_to_cnsq(int &cnsq,
		std::vector<int> const& xyz,
		std::vector<int> const& grid,
		bool dataLayoutRowMajor)
{
	assert( xyz.size() == grid.size() );
	cnsq = 0;
	int dim = grid.size();
	if ( dataLayoutRowMajor )
	{
		//e.g. for d == 3 : ix*d[1]+iy)*d[2]+iz
		for ( int i = 0 ; i < dim-1; ++i)
		{
			cnsq += xyz[i];
			cnsq *= grid[i+1];
		}
		int im = dim-1;
		cnsq += xyz[im];
	}
	else // Column major
	{
		//e.g. for d == 3 : ix+d[0]*(iy+d[1]*iz
		for ( int i = dim-1 ; i > 0; --i)
		{
			cnsq += xyz[i];
			cnsq *= grid[i-1];
		}
		cnsq += xyz[0];
	}
}

} /* namespace Algorithms */
} /* namespace elephon */
