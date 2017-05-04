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
		std::vector<size_t> grid,
		std::vector<double> const& latticeMatrix,
		size_t nDataBlock,
		std::vector<double> const& dataOnGrid)
{
	grid_ = std::move( grid );
	nBlockData_ = nDataBlock;
	size_t gridnum = grid_[0]*grid_[1]*grid_[2];
	assert(latticeMatrix.size()==9);
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
	for ( int id = 0 ; id < dim; ++id)
		n[id] = grid_[id];
	int dblock = gridnum;
	int * inembed =NULL;
	int * onembed = NULL;
	int howmany = static_cast<int>( nBlockData_ );
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
	howmany = static_cast<int>( nBlockData_*3 );

	fftw_complex * FFTBufferRtoK = fftw_alloc_complex( dataOnGrid.size()*3 );
	auto fftw3PlanFwdRtoK = fftw_plan_many_dft(
			dim,n,howmany,
			FFTBufferRtoK, inembed,istride, idist,
			FFTBufferRtoK, onembed,ostride, odist,
			-1,FFTW_ESTIMATE);

	for ( size_t ig = 0 ; ig < gridnum; ig++)
		for ( size_t id = 0 ; id < nBlockData_; id++)
			reinterpret_cast<std::complex<double>*>(FFTBuffer)[id*gridnum+ig]
					= std::complex<double>(dataOnGrid[ig*nBlockData_+id]);

	fftw_execute( fftw3PlanBkwdKtoR );

	//define a lambda to transform an consecutively ordered index into a Cartesian R vector
	//All grids are z fastest running and x slowest
	std::vector<int> indexTouple(3);

	//For the derivative in direction x we need to have the strictly balanced FT
	//Thus if the grid is even, we need to set the unbalanced largest component to zero.
	bool setZero[3] = { not static_cast<bool>(grid_[0]%2),
			not static_cast<bool>(grid_[1]%2),
			not static_cast<bool>(grid_[2]%2) };

	std::vector<double> RCartesianVectors(3*gridnum);
	for (size_t ig = 0 ; ig < gridnum; ++ig)
	{
		//reverse-engineer expressions
		// index = (x*dy + y)*dz + z
		//		  = x*dy*dz + y*dz + z
		//for x,y,z using integer division

		int d = 3;
		std::vector<int> dimProd( d , 1 );
		for (int i = d-2 ; i >= 0 ; --i)
			dimProd[i] = grid_[ i+1 ]*dimProd[i+1];

	 	indexTouple.assign( d , ig );
		for ( int i = 0 ;  i < d-1 ; ++i )
		{
			indexTouple[i] = indexTouple[i]/dimProd[i];
			for ( int j = i+1 ;  j < d ; ++j )
				indexTouple[j] -= indexTouple[i]*dimProd[i];
		}

		//CAREFUL! We need to perform this strictly in the 1st unit cell!
		for ( int i = 0 ;  i < 3 ; ++i )
		{
			if ( setZero[i] and (indexTouple[i] == static_cast<int>(grid_[i]/2)) )
				indexTouple[i] = 0;

			if ( indexTouple[i] > static_cast<int>(grid_[i]/2) )
				indexTouple[i] -= static_cast<int>(grid_[i]);
		}

		//Transform the grid vector into Cartesian coordinates by
		//multiplying the lattice matrix
		for ( size_t i = 0 ; i < 3 ; i++)
		{
			RCartesianVectors[ig*3+i] = 0 ;
			for ( size_t j = 0 ; j < 3 ; j++)
			{
				double xj = double(indexTouple[j]);
				RCartesianVectors[ig*3+i] += xj*latticeMatrix[i*3+j];
			}
		}
	};

	for ( size_t ig = 0 ; ig < gridnum; ig++)
		for ( size_t id = 0 ; id < nBlockData_; id++)
			for (size_t ix = 0 ; ix < 3 ; ix++)
				reinterpret_cast<std::complex<double>*>(FFTBufferRtoK)[(id*3+ix)*gridnum+ig]
                  = -2.0*M_PI*std::complex<double>(0,RCartesianVectors[ig*3+ix])*
				  reinterpret_cast<std::complex<double>*>(FFTBuffer)[id*gridnum+ig];

	fftw_execute( fftw3PlanFwdRtoK );

	gradientDataOnGrid_ = std::vector< double >( gridnum*nBlockData_*3 );

	double imagCheck = 0;
	for ( size_t ig = 0 ; ig < gridnum; ig++)
		for ( size_t id = 0 ; id < nBlockData_; id++)
			for ( size_t ix = 0 ; ix < 3; ix++)
			{
				gradientDataOnGrid_[(ig*nBlockData_+id)*3+ix] =
						std::real(reinterpret_cast<std::complex<double>*>(
								FFTBufferRtoK)[(id*3+ix)*gridnum+ig])
						/double(gridnum);//normalization
				imagCheck +=
						std::imag(reinterpret_cast<std::complex<double>*>(
								FFTBufferRtoK)[(id*3+ix)*gridnum+ig]);
			}
	assert(imagCheck/double(gridnum) < 1e-6);

	delete [] n;
	delete [] inembed;
	delete [] onembed;

	fftw_free( FFTBuffer );
	fftw_free( FFTBufferRtoK );
}

void GradientFFTReciprocalGrid::copy_data(std::vector<size_t> const& conseqGridIndices,
		std::vector<size_t> bands,
		std::vector<double> & dataAtIndices) const
{
	if ( dataAtIndices.size() != conseqGridIndices.size()*bands.size()*3 )
		dataAtIndices = std::vector<double>(conseqGridIndices.size()*bands.size()*3 );
	for (size_t i = 0 ; i < conseqGridIndices.size(); ++i )
		for ( size_t id = 0 ; id < bands.size(); id++)
		{
			size_t iband = bands[id];
			assert( iband < nBlockData_ );
			for ( size_t ix = 0 ; ix < 3; ix++)
				dataAtIndices[(i*bands.size()+id)*3+ix] = gradientDataOnGrid_[ (conseqGridIndices[i]*nBlockData_+iband)*3+ix ];
		}
}

std::vector<double> const& GradientFFTReciprocalGrid::get_data() const
{
	return gradientDataOnGrid_;
}

} /* namespace ElectronicStructure */
} /* namespace elephon */
