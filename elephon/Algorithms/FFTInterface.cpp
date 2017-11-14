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
#include <omp.h>

namespace elephon
{
namespace Algorithms
{

FFTInterface::FFTInterface()
{

}

void
FFTInterface::plan_fft_local(
		std::vector<int> const & gridDims,
		int exponentSign,
		int nDataPerGridPt,
		int hintHowOften,
		fftw_plan_s * & toAlloc,
		bool dataLayoutRowMajor)
{
	if ( omp_get_num_threads() != 1 )
		throw std::runtime_error("Called FFTInterface::plan_fft in a parallel region: not thread safe");

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

	// please note: the planner takes always the first FFTBuffer. We call with different fftw_exectute
	// as described at http://www.fftw.org/fftw3_doc/New_002darray-Execute-Functions.html#New_002darray-Execute-Functions
	toAlloc = fftw_plan_many_dft(
				dim,n,howmany,
				FFTBuffer_[0], inembed,istride, idist,
				FFTBuffer_[0], onembed,ostride, odist,
				exponentSign,plan_how_many);

	delete [] n;
	delete [] inembed;
	delete [] onembed;
}

FFTInterface::~FFTInterface()
{
	this->clear_storadge();
}

void
FFTInterface::clear_storadge()
{
	dataLayoutRowMajor_ = false;
	nDataPerGridPt_ = 1;
	gridDimsBKWD_.clear();
	gridDimsFWD_.clear();
	nBuff_.clear();
	if ( FFTBuffer_ != nullptr )
	{
		for (int it = 0 ; it < numThreadsMax_; ++it)
			fftw_free( FFTBuffer_[it] );
		delete [] FFTBuffer_;
		FFTBuffer_ = nullptr;
	}
	if ( fftw3PlanFowdRtoK_ != nullptr )
	{
		fftw_destroy_plan( fftw3PlanFowdRtoK_ );
		fftw3PlanFowdRtoK_ = nullptr;
	}
	if ( fftw3PlanBkwdKtoR_ != nullptr )
	{
		fftw_destroy_plan( fftw3PlanBkwdKtoR_ );
		fftw3PlanBkwdKtoR_ = nullptr;
	}
}

void
FFTInterface::plan_fft(
		std::vector<int> const & gridDimsData,
		int nDataPerGridPt,
		int exponentSign,
		bool dataLayoutRowMajor,
		int hintHowOften)
{
	if ( omp_get_num_threads() != 1 )
		throw std::runtime_error("Error in FFTInterface::plan_fft: not thread safe but called with #threads != 1");
	this->allocate_internal();

	dataLayoutRowMajor_ = dataLayoutRowMajor;
	nDataPerGridPt_ = nDataPerGridPt;

	this->allocate_this_thread(gridDimsData, nDataPerGridPt_); // will allocate thread 0

	if ( ! (exponentSign < 0) ) //  > 0 and 0
	{
		if (gridDimsData != gridDimsBKWD_)
		{
			gridDimsBKWD_ = gridDimsData;
			this->plan_fft_local(gridDimsBKWD_, +1, nDataPerGridPt_, hintHowOften, fftw3PlanBkwdKtoR_, dataLayoutRowMajor);
		}
	}
	if (  ! (exponentSign > 0) ) //  < 0 and 0
	{
		if (gridDimsData != gridDimsFWD_)
		{
			gridDimsFWD_ = gridDimsData;
			this->plan_fft_local(gridDimsFWD_, -1, nDataPerGridPt_, hintHowOften, fftw3PlanFowdRtoK_ , dataLayoutRowMajor);
		}
	}
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

void
FFTInterface::allocate_this_thread(
		std::vector<int> const & gridDims,
		int nDataPerGridPt)
{
	if (  FFTBuffer_ == nullptr  )
		throw std::logic_error("FFTInterface::allocate_this_thread: must call allocate_internal before allocate this thread");
	int threadID = omp_get_thread_num();
	assert(threadID < numThreadsMax_);

	int ngrid = 1;
	for ( auto di : gridDims )
		ngrid *= di;
	int nTotal = nDataPerGridPt*ngrid;

	assert(nBuff_.size() == numThreadsMax_);
	if ( nTotal > nBuff_[threadID] )
	{
		if (FFTBuffer_[threadID] != nullptr)
		{
			fftw_free( FFTBuffer_[threadID] );
			FFTBuffer_[threadID] = nullptr;
		}
		FFTBuffer_[threadID] = fftw_alloc_complex( nTotal );
		nBuff_[threadID] = nTotal;
	}
}

void
FFTInterface::allocate_internal()
{
	numThreadsMax_ = omp_get_max_threads();
	if ( omp_get_num_threads() != 1 )
		throw std::runtime_error("Error in FFTInterface::allocate_internal: not thread safe but called with #threads != 1");

	this->clear_storadge();

	if ( FFTBuffer_ == nullptr )
	{
		FFTBuffer_ = new fftw_complex* [numThreadsMax_];
		for (int it = 0 ; it < numThreadsMax_; ++it)
			FFTBuffer_[it] = nullptr;
	}

	nBuff_.assign(numThreadsMax_, 0);
}

} /* namespace Algorithms */
} /* namespace elephon */
