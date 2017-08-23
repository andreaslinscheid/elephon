/*	This file FFTInterface.hpp is part of elephon.
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
#include <stdexcept>
#include <iostream>

namespace elephon
{
namespace Algorithms
{

namespace detail
{

template<class To, class Ti>
struct ComplexConversion
{
	To convert(Ti val) {
		return To(val);
	}

	const bool complex_to_real = false;

	//Not used here
	To imagAccumalate = To(0);
};

template<class To, class Ti>
struct ComplexConversion<std::complex<To>,std::complex<Ti> >
{
	std::complex<To> convert(std::complex<Ti> val) {
		return std::complex<To>(val);
	}

	const bool complex_to_real = false;

	//Not used here
	To imagAccumalate = To(0);
};

template<class To, class Ti>
struct ComplexConversion<To,std::complex<Ti> >
{
	To convert(std::complex<Ti> val) {
		imagAccumalate += std::abs(std::imag(val));
		return To(std::real(val));
	}

	const bool complex_to_real = true;

	To imagAccumalate = To(0);
};

template<class T>
struct MakeComplex
{
	typedef std::complex<T> type;
};

template<class T>
struct MakeComplex<std::complex<T>>
{
	typedef std::complex<T> type;
};
} /* namespace detail */

template< typename TR>
void
FFTInterface::allocate(
		std::vector<int> const & gridDims,
		int nDataPerGridPt,
		std::vector< TR > & dataResult,
		bool & re_plan)
{
	int ngrid = 1;
	for ( auto di : gridDims )
		ngrid *= di;
	int nTotal = nDataPerGridPt*ngrid;

	re_plan = false;

	if ( nTotal > nBuff_ )
	{
		if (FFTBuffer_ != 0)
		{
			fftw_free( FFTBuffer_ );
			FFTBuffer_ = nullptr;
		}
		FFTBuffer_ = fftw_alloc_complex( nTotal );
		nBuff_ = nTotal;
		re_plan = true;
	}
	dataResult.resize( nTotal );
}

template<typename TI, typename TR>
void
FFTInterface::fft_sparse_data(
		std::vector<int> const & mapFFTCoeff,
		std::vector<int> const & gridDimsInputData,
		std::vector< TI > const & sparseInputData,
		int nDataPerGridPt,
		int exponentSign,
		std::vector< TR > & dataResult,
		std::vector<int> const & gridDimsOutputData,
		bool dataLayoutRowMajor,
		int hintHowOften )
{
	int spaceDim = gridDimsOutputData.size();
	assert(mapFFTCoeff.size()/spaceDim == sparseInputData.size()/nDataPerGridPt);
	int ngrid = 1;
	for ( auto di : gridDimsOutputData )
		ngrid *= di;

	bool re_plan = false;;
	this->allocate(gridDimsOutputData, nDataPerGridPt, dataResult, re_plan);

	detail::ComplexConversion< std::complex<double>, TI > converterFwd;

	auto xyz_to_cnsq = [&] (
			std::vector<int>::const_iterator toupleBegin,
			std::vector<int>::const_iterator toupleEnd)
	{
		auto output_aligned = [] (int iG, int maxGIn, int maxGOut)
		{
			//Convert to the negative-positive frequency scheme
			int iGAct = iG < maxGIn/2 ? iG : iG - maxGIn;
			//alias back in case ...
			iGAct = iGAct % maxGOut;
			//Convert back to the positive index scheme with the second half being the negative freqs
			return iGAct < 0 ? iGAct + maxGOut : iGAct;
		};
		assert( std::distance(toupleBegin,toupleEnd) == gridDimsOutputData.size() );
		auto const & d = gridDimsOutputData;
		int conseq = 0;
		if ( dataLayoutRowMajor )
		{
			//e.g. for d == 3 : ix*d[1]+iy)*d[2]+iz
			for ( int i = 0 ; i < int(d.size())-1; ++i)
			{
				conseq += output_aligned(*(toupleBegin+i), gridDimsInputData[i], d[i]);
				conseq *= d[i+1];
			}
			int im = int(d.size())-1;
			conseq += output_aligned(*(toupleBegin+im), gridDimsInputData[im], d[im]);
		}
		else // Column major
		{
			//e.g. for d == 3 : ix+d[0]*(iy+d[1]*iz
			for ( int i = int(d.size())-1 ; i > 0; --i)
			{
				conseq += output_aligned(*(toupleBegin+i), gridDimsInputData[i], d[i]);
				conseq *= d[i-1];
			}
			conseq += output_aligned(*(toupleBegin), gridDimsInputData[0], d[0]);
		}
		return conseq;
	};

	auto fillBuffer = [&] () {
		//Not all data may be accessed, so set it to zero explictely
		auto ptr = reinterpret_cast<std::complex<double> * >(FFTBuffer_);
		std::fill(ptr,ptr+ngrid*nDataPerGridPt, std::complex<double>(0));

		int nGridPointsNonZero = mapFFTCoeff.size()/spaceDim;
		for (int ig = 0 ; ig < nGridPointsNonZero; ++ig)
		{
			auto it = mapFFTCoeff.begin()+ig*spaceDim;
			int conseq_index = xyz_to_cnsq(it,it+spaceDim);
			for (int id = 0 ; id < nDataPerGridPt; ++id)
			{
				ptr[ id*ngrid + conseq_index ] = converterFwd.convert( sparseInputData [id*nGridPointsNonZero + ig] );
				assert( ptr[ id*ngrid + conseq_index ] == ptr[ id*ngrid + conseq_index ] );
			}
		}
	};

	auto fillResult = [&] () {
		detail::ComplexConversion< TR, std::complex<double> > converterBkwd;
		for (int id = 0 ; id < nDataPerGridPt; ++id)
			for (int ig = 0 ; ig < ngrid; ++ig)
			{
				dataResult [ id*ngrid + ig] =
						converterBkwd.convert( reinterpret_cast<std::complex<double> * >(FFTBuffer_)[ id*ngrid + ig ] );
				assert( dataResult [ id*ngrid + ig] == dataResult [ id*ngrid + ig] );
			}

		if ( converterBkwd.complex_to_real )
			if ( std::abs(converterBkwd.imagAccumalate)/ngrid/nDataPerGridPt > 1e-6 )
				throw std::runtime_error(
						std::string("FFT from complex to real lost significant information by slicing an imaginary part of ")
						+std::to_string(std::abs(converterBkwd.imagAccumalate)/ngrid/nDataPerGridPt) +
						" on average per data point");

	};

	//the grid dimension are the flag to recompute a plan, both in forward and backward direction
	if ( exponentSign > 0 )
	{
		if ( (gridDimsOutputData != gridDimsBKWD_) or re_plan )
		{
			gridDimsBKWD_ = gridDimsOutputData;
			this->plan_fft(gridDimsBKWD_, +1, nDataPerGridPt, hintHowOften, fftw3PlanBkwdKtoR_, dataLayoutRowMajor);
		}

		fillBuffer();
		fftw_execute(fftw3PlanBkwdKtoR_);
		fillResult();
	}
	else
	{
		if ( (gridDimsOutputData != gridDimsFWD_) or re_plan )
		{
			gridDimsFWD_ = gridDimsOutputData;
			this->plan_fft(gridDimsFWD_, -1, nDataPerGridPt, hintHowOften, fftw3PlanFowdRtoK_ , dataLayoutRowMajor);
		}

		fillBuffer();
		fftw_execute(fftw3PlanFowdRtoK_);
		fillResult();
	}
};


template< typename TR>
void
FFTInterface::fill_result(int ngrid, int nDataPerGridPt, std::vector<TR> & dataResult )
{
	detail::ComplexConversion< TR, std::complex<double> > converterBkwd;
	for (int ig = 0 ; ig < ngrid; ++ig)
		for (int id = 0 ; id < nDataPerGridPt; ++id)
			dataResult [ ig * nDataPerGridPt + id ] =
					converterBkwd.convert( reinterpret_cast<std::complex<double> * >(FFTBuffer_)[ id*ngrid + ig ] );
	if ( converterBkwd.complex_to_real )
		if ( std::abs(converterBkwd.imagAccumalate)/ngrid/nDataPerGridPt > 1e-6 )
			throw std::runtime_error(
					std::string("FFT from complex to real lost significant information by slicing an imaginary part of ")
					+std::to_string(std::abs(converterBkwd.imagAccumalate)/ngrid/nDataPerGridPt) + " on average per data point");
};

template<typename TI, typename TR>
void
FFTInterface::fft_data(
		std::vector<int> const & gridDims,
		std::vector< TI > const & data,
		std::vector< TR > & dataResult,
		int nDataPerGridPt,
		int exponentSign,
		bool dataLayoutRowMajor,
		int hintHowOften )
{
	int ngrid = 1;
	for ( auto di : gridDims )
		ngrid *= di;

	bool re_plan = false;;
	this->allocate(gridDims, nDataPerGridPt, dataResult, re_plan);

	detail::ComplexConversion< std::complex<double>, TI > converterFwd;

	auto fillBuffer = [&] () {
		for (int ig = 0 ; ig < ngrid; ++ig)
			for (int id = 0 ; id < nDataPerGridPt; ++id)
				reinterpret_cast<std::complex<double> * >(FFTBuffer_)[ id*ngrid + ig ] =
						converterFwd.convert( data [ ig * nDataPerGridPt + id ] );
	};

	//the grid dimension are the flag to recompute a plan, both in forward and backward direction
	if ( exponentSign > 0 )
	{
		if ( (gridDims != gridDimsBKWD_) or re_plan )
		{
			gridDimsBKWD_ = gridDims;
			this->plan_fft(gridDimsBKWD_, +1, nDataPerGridPt, hintHowOften, fftw3PlanBkwdKtoR_, dataLayoutRowMajor);
		}

		fillBuffer();
		fftw_execute(fftw3PlanBkwdKtoR_);
		this->fill_result(ngrid,nDataPerGridPt,dataResult);
	}
	else
	{
		if ( (gridDims != gridDimsFWD_) or re_plan )
		{
			gridDimsFWD_ = gridDims;
			this->plan_fft(gridDimsFWD_, -1, nDataPerGridPt, hintHowOften, fftw3PlanFowdRtoK_ , dataLayoutRowMajor);
		}

		fillBuffer();
		fftw_execute(fftw3PlanFowdRtoK_);
		this->fill_result(ngrid,nDataPerGridPt,dataResult);
	}
}

template<typename T>
void
FFTInterface::fft_interpolate(
               std::vector<int> const & gridDimsIn,
               std::vector< T > const & data,
               std::vector<int> const & gridDimsOut,
               std::vector< T > & dataResult,
               int nDataPerGridPt)
{
       assert(gridDimsIn.size() == gridDimsOut.size());
       assert(gridDimsIn.size() == 3);
       int ngridIn = gridDimsIn[0]*gridDimsIn[1]*gridDimsIn[2];
       int ngridOut = gridDimsOut[0]*gridDimsOut[1]*gridDimsOut[2];
       assert( (ngridIn != 0) && (ngridOut != 0));

       // here we assume x major grid order
       std::vector< typename detail::MakeComplex<T>::type > intermedOut;
       this->fft_data(gridDimsIn, data, intermedOut, nDataPerGridPt, -1, false, 1);

       for ( auto &d : intermedOut )
               d /= T(ngridIn);

       std::vector< typename detail::MakeComplex<T>::type > intermedIn(ngridOut*nDataPerGridPt,
                                                               typename detail::MakeComplex<T>::type(0.0));
       for ( int iz = 0 ; iz < gridDimsIn[2]; ++iz)
               for ( int iy = 0 ; iy < gridDimsIn[1]; ++iy)
                       for ( int ix = 0 ; ix < gridDimsIn[0]; ++ix)
                       {
                               int cnsq = ix + gridDimsIn[0]*(iy + gridDimsIn[1]*iz);
                               // go in to the frequency representation with + and - frequencies in the input mesh
                               int iX = ix;
                               int iY = iy;
                               int iZ = iz;
                               Algorithms::FFTInterface::inplace_to_freq(iX, iY, iZ, gridDimsIn[0], gridDimsIn[1], gridDimsIn[2]);
                               // frequencies outside the output grid are dropped
                               if ( ((iX < -gridDimsOut[0]/2-gridDimsOut[0]%2) or (iX >= gridDimsOut[0]/2+gridDimsOut[0]%2)) or
                                        ((iY < -gridDimsOut[1]/2-gridDimsOut[1]%2) or (iY >= gridDimsOut[1]/2+gridDimsOut[1]%2)) or
                                        ((iZ < -gridDimsOut[2]/2-gridDimsOut[2]%2) or (iZ >= gridDimsOut[2]/2+gridDimsOut[2]%2)) )
                                       continue;
                               // go in to the fftw3 representation in the positive frequency notation in the output mesh
                               Algorithms::FFTInterface::freq_to_inplace(iX, iY, iZ, gridDimsOut[0], gridDimsOut[1], gridDimsOut[2]);
                               int cnsqNew = iX + gridDimsOut[0]*(iY + gridDimsOut[1]*iZ);
                               for ( int ib = 0 ; ib < nDataPerGridPt; ++ib )
                                       intermedIn[cnsqNew*nDataPerGridPt+ib] = intermedOut[cnsq*nDataPerGridPt+ib];
                       }

       this->fft_data(gridDimsOut, intermedIn, dataResult, nDataPerGridPt, 1, false, 1);
}


} /* namespace Algorithms */
} /* namespace elephon */
