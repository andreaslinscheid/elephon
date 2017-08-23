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
	FFTInterface::inplace_to_freq(g[0], g[1], g[2], d[0], d[1], d[2]);
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
	FFTInterface::freq_to_inplace(g[0], g[1], g[2], d[0], d[1], d[2]);
}

} /* namespace Algorithms */
} /* namespace elephon */
