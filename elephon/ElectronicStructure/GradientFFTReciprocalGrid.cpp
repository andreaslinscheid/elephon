/*	This file GradientFFTReciprocalGrid.cpp is part of elephon.
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
 *  Created on: Apr 29, 2017
 *      Author: A. Linscheid
 */

#include "GradientFFTReciprocalGrid.h"
#include <fftw3.h>
#include <complex>
#include <assert.h>
#include <iostream>

namespace elephon
{
namespace ElectronicStructure
{

GradientFFTReciprocalGrid::GradientFFTReciprocalGrid()
{

}

void GradientFFTReciprocalGrid::compute_gradient(
		std::vector<int> grid,
		LatticeStructure::LatticeModule const& lattice,
		int nDataBlock,
		std::vector<double> const& dataOnGrid)
{
	grid_ = std::move( grid );
	nBlockData_ = nDataBlock;
	int gridnum = grid_[0]*grid_[1]*grid_[2];
	assert(dataOnGrid.size() == gridnum*nBlockData_);

	//Plan the FFT:
	//We provide fftw3 the data for many DFTs laid out as each one band (and all 3 directions) for the grid.
	//This means we need to shuffle the data around ...
	fftw_complex * FFTBuffer = fftw_alloc_complex( dataOnGrid.size() );

	//Please see the fftw3 documentation for the meaning of the variables below.
	//I am using the same nomenclature.
	//TODO This can be done slightly better using the real to complex FFT
	int dim = 3;
	int * n = new int [ dim ];
	n[0] = grid_[0];
	n[1] = grid_[1];
	n[2] = grid_[2];
	int dblock = gridnum;
	int * inembed =NULL;
	int * onembed = NULL;
	int howmany = nBlockData_;
	int idist = dblock;
	int odist = dblock;
	int istride = 1;
	int ostride = 1;

	auto fftw3PlanBkwdKtoR = fftw_plan_many_dft(
			dim,n,howmany,
			FFTBuffer, inembed,istride, idist,
			FFTBuffer, onembed,ostride, odist,
			1,FFTW_ESTIMATE);

	//We multiply i R in Cartesian coordinates and transform back to obtain the gradient
	howmany = nBlockData_*3 ;

	fftw_complex * FFTBufferRtoK = fftw_alloc_complex( dataOnGrid.size()*3 );
	auto fftw3PlanFwdRtoK = fftw_plan_many_dft(
			dim,n,howmany,
			FFTBufferRtoK, inembed,istride, idist,
			FFTBufferRtoK, onembed,ostride, odist,
			-1,FFTW_ESTIMATE);

	for ( int iGz = 0 ; iGz < grid_[2]; ++iGz )
		for ( int iGy = 0 ; iGy < grid_[1]; ++iGy )
			for ( int iGx = 0 ; iGx < grid_[0]; ++iGx )
			{
				int fftw3Layout = iGz + grid_[2]*(iGy+grid_[1]*iGx);
				int presentLayout = iGx + grid_[0]*(iGy+grid_[1]*iGz);
				for ( int id = 0 ; id < nBlockData_; id++)
				{
					reinterpret_cast<std::complex<double>*>(FFTBuffer)[id*gridnum+fftw3Layout]
							= std::complex<double>(dataOnGrid[presentLayout*nBlockData_+id]);
				}
			}

	fftw_execute( fftw3PlanBkwdKtoR );

	//All grids are x fastest running and z slowest
	//For the derivative in direction x we need to have the strictly balanced FT
	//Thus if the grid is even, we need to set the unbalanced largest component to zero.
	bool setZero[3] = { not static_cast<bool>(grid_[0]%2),
			not static_cast<bool>(grid_[1]%2),
			not static_cast<bool>(grid_[2]%2) };

	//Actually, we are not setting R to zero, but we are simply not assigning anything else than
	//zero in the code below. This multiplies the data, effectively dropping this contribution.
	std::vector<double> RCartesianVectors( 3*gridnum , 0.0 );
	for ( int iGz = 0 ; iGz < grid_[2]; ++iGz )
	{
		int iGzf = iGz < grid_[2]/2 ? iGz : iGz-grid_[2];
		if ( setZero[2] and ( iGz == grid_[2]/2) )
			continue;
		for ( int iGy = 0 ; iGy < grid_[1]; ++iGy )
		{
			int iGyf = iGy < grid_[1]/2 ? iGy : iGy-grid_[1];
			if ( setZero[1] and ( iGy == grid_[1]/2) )
				continue;
			for ( int iGx = 0 ; iGx < grid_[0]; ++iGx )
			{
				int iGxf = iGx < grid_[0]/2 ? iGx : iGx-grid_[0];
				if ( setZero[0] and ( iGx == grid_[0]/2) )
					continue;

				//NOTE: R multiples fftw3 laid out data which is z fastest running!
				int ig = (iGx*grid_[1]+iGy)*grid_[2]+iGz;

				RCartesianVectors[ig*3+0] = iGxf;
				RCartesianVectors[ig*3+1] = iGyf;
				RCartesianVectors[ig*3+2] = iGzf;
			}
		}
	}

	lattice.direct_to_cartesian_angstroem(RCartesianVectors);

	for ( int ig = 0 ; ig < gridnum; ig++)
		for ( int id = 0 ; id < nBlockData_; id++)
			for (int ix = 0 ; ix < 3 ; ix++)
				reinterpret_cast<std::complex<double>*>(FFTBufferRtoK)[(id*3+ix)*gridnum+ig]
                  = -std::complex<double>(0,RCartesianVectors[ig*3+ix])*
				  reinterpret_cast<std::complex<double>*>(FFTBuffer)[id*gridnum+ig];

	fftw_execute( fftw3PlanFwdRtoK );

	gradientDataOnGrid_ = std::vector< double >( gridnum*nBlockData_*3 );

	double imagCheck = 0;	for ( int iGz = 0 ; iGz < grid_[2]; ++iGz )
	for ( int iGy = 0 ; iGy < grid_[1]; ++iGy )
		for ( int iGx = 0 ; iGx < grid_[0]; ++iGx )
		{
			int fftw3Layout = iGz + grid_[2]*(iGy+grid_[1]*iGx);
			int presentLayout = iGx + grid_[0]*(iGy+grid_[1]*iGz);
			for ( int id = 0 ; id < nBlockData_; id++)
				for ( int ix = 0 ; ix < 3; ix++)
				{
					gradientDataOnGrid_[(presentLayout*nBlockData_+id)*3+ix] =
							std::real(reinterpret_cast<std::complex<double>*>(
									FFTBufferRtoK)[(id*3+ix)*gridnum+fftw3Layout])
							/double(gridnum);//normalization
					imagCheck +=
							std::imag(reinterpret_cast<std::complex<double>*>(
									FFTBufferRtoK)[(id*3+ix)*gridnum+fftw3Layout]);
				}
		}
	assert(imagCheck/double(gridnum) < 1e-6);

	delete [] n;
	delete [] inembed;
	delete [] onembed;

	fftw_free( FFTBuffer );
	fftw_free( FFTBufferRtoK );
}

void GradientFFTReciprocalGrid::copy_data(std::vector<int> const& conseqGridIndices,
		std::vector<int> bands,
		std::vector<double> & dataAtIndices) const
{
	if ( dataAtIndices.size() != conseqGridIndices.size()*bands.size()*3 )
		dataAtIndices.resize(conseqGridIndices.size()*bands.size()*3 );
	for (int i = 0 ; i < conseqGridIndices.size(); ++i )
		for ( int id = 0 ; id < bands.size(); id++)
		{
			int iband = bands[id];
			assert( iband < nBlockData_ );
			for ( int ix = 0 ; ix < 3; ix++)
				dataAtIndices[(i*bands.size()+id)*3+ix] =
						gradientDataOnGrid_[ (conseqGridIndices[i]*nBlockData_+iband)*3+ix ];
		}
}

std::vector<double> const& GradientFFTReciprocalGrid::get_data() const
{
	return gradientDataOnGrid_;
}

} /* namespace ElectronicStructure */
} /* namespace elephon */
