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
#include <algorithm>
#include <omp.h>

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

template<typename VI, typename VR>
void
FFTInterface::fft_sparse_data(
		std::vector<int> const & mapFFTCoeff,
		std::vector<int> const & gridDimsInputData,
		VI const & sparseInputData,
		int exponentSign,
		VR & dataResult )
{
	assert(exponentSign != 0);
	int threadID = omp_get_thread_num();
	if ( threadID >= numThreadsMax_)
		throw std::runtime_error("FFTInterface::fft_sparse_data: the number of threads is beyond what it was when initialized");

	std::vector<int> gridDimsOutputData = exponentSign > 0 ? gridDimsBKWD_ : gridDimsFWD_;
	if (gridDimsOutputData.size() == 0)
		throw std::runtime_error("FFTInterface::fft_sparse_data need to plan direction first");

	int spaceDim = gridDimsOutputData.size();
	assert(mapFFTCoeff.size()/spaceDim == sparseInputData.size()/nDataPerGridPt_);
	int ngrid = 1;
	for ( auto di : gridDimsOutputData )
		ngrid *= di;

	this->allocate_this_thread(gridDimsOutputData, nDataPerGridPt_);

	dataResult.resize(ngrid*nDataPerGridPt_);

	detail::ComplexConversion< std::complex<double>, typename VI::value_type > converterFwd;

	auto xyz_to_cnsq = [&] (
			std::vector<int>::const_iterator toupleBegin,
			std::vector<int>::const_iterator toupleEnd)
	{
		auto output_aligned = [] (int iG, int maxGIn, int maxGOut)
		{
			//Convert to the negative-positive frequency scheme
			int iGAct = iG <= maxGIn/2 ? iG : iG - maxGIn;
			//alias back in case ...
			iGAct = iGAct % maxGOut;
			//Convert back to the positive index scheme with the second half being the negative freqs
			return iGAct < 0 ? iGAct + maxGOut : iGAct;
		};
		assert( std::distance(toupleBegin,toupleEnd) == gridDimsOutputData.size() );
		auto const & d = gridDimsOutputData;
		int conseq = 0;
		if ( dataLayoutRowMajor_ )
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
		assert(FFTBuffer_ != nullptr);

		//Not all data may be accessed, so set it to zero explictely
		auto ptr = reinterpret_cast<std::complex<double> * >(FFTBuffer_[threadID]);
		std::fill(ptr,ptr+ngrid*nDataPerGridPt_, std::complex<double>(0));

		int nGridPointsNonZero = mapFFTCoeff.size()/spaceDim;
		for (int ig = 0 ; ig < nGridPointsNonZero; ++ig)
		{
			auto it = mapFFTCoeff.begin()+ig*spaceDim;
			int conseq_index = xyz_to_cnsq(it,it+spaceDim);
			for (int id = 0 ; id < nDataPerGridPt_; ++id)
			{
				ptr[ id*ngrid + conseq_index ] = converterFwd.convert( sparseInputData [id*nGridPointsNonZero + ig] );
				assert( ptr[ id*ngrid + conseq_index ] == ptr[ id*ngrid + conseq_index ] );
			}
		}
	};

	auto fillResult = [&] () {
		assert(FFTBuffer_ != nullptr);

		detail::ComplexConversion< typename VR::value_type, std::complex<double> > converterBkwd;
		for (int id = 0 ; id < nDataPerGridPt_; ++id)
			for (int ig = 0 ; ig < ngrid; ++ig)
			{
				dataResult [ id*ngrid + ig] =
						converterBkwd.convert( reinterpret_cast<std::complex<double> * >(FFTBuffer_[threadID])[ id*ngrid + ig ] );
				assert( dataResult [ id*ngrid + ig] == dataResult [ id*ngrid + ig] );
			}

		if ( converterBkwd.complex_to_real )
			if ( std::abs(converterBkwd.imagAccumalate)/ngrid/nDataPerGridPt_ > 1e-6 )
			{
				auto realIne = decltype(converterBkwd.imagAccumalate)(0);
				for (int id = 0 ; id < nDataPerGridPt_; ++id)
					for (int ig = 0 ; ig < ngrid; ++ig)
						realIne += std::real(dataResult [ id*ngrid + ig]);
				if ( std::abs(converterBkwd.imagAccumalate)/std::abs(realIne) > 1e-6 )
					throw std::runtime_error(
							std::string("FFT from complex to real lost significant information by slicing an imaginary part of ")
							+std::to_string(std::abs(converterBkwd.imagAccumalate)/ngrid/nDataPerGridPt_) +
							" on average per data point");
			}
	};

	if ( exponentSign > 0 )
	{
		fillBuffer();
		fftw_execute_dft(fftw3PlanBkwdKtoR_, FFTBuffer_[threadID], FFTBuffer_[threadID]);
		fillResult();
	}
	else
	{
		fillBuffer();
		fftw_execute_dft(fftw3PlanFowdRtoK_, FFTBuffer_[threadID], FFTBuffer_[threadID]);
		fillResult();
	}
};


template< typename VTR>
void
FFTInterface::fill_result(int ngrid, VTR & dataResult )
{
	int threadID = omp_get_thread_num();
	dataResult.resize(ngrid*nDataPerGridPt_);
	detail::ComplexConversion< typename VTR::value_type, std::complex<double> > converterBkwd;
	for (int ig = 0 ; ig < ngrid; ++ig)
		for (int id = 0 ; id < nDataPerGridPt_; ++id)
		{
			dataResult [ ig * nDataPerGridPt_ + id ] =
					converterBkwd.convert( reinterpret_cast<std::complex<double> * >(FFTBuffer_[threadID])[ id*ngrid + ig ] );
			assert( dataResult [ ig * nDataPerGridPt_ + id ] == dataResult [ ig * nDataPerGridPt_ + id ]);
		}
	if ( converterBkwd.complex_to_real )
		if ( std::abs(converterBkwd.imagAccumalate)/ngrid/nDataPerGridPt_ > 1e-6 )
		{
			auto realIne = decltype(converterBkwd.imagAccumalate)(0);
			for (int id = 0 ; id < nDataPerGridPt_; ++id)
				for (int ig = 0 ; ig < ngrid; ++ig)
					realIne += std::real(dataResult [ id*ngrid + ig]);
			if ( std::abs(converterBkwd.imagAccumalate)/std::abs(realIne) > 1e-6 )
				throw std::runtime_error(
						std::string("FFT from complex to real lost significant information by slicing an imaginary part of ")
						+std::to_string(std::abs(converterBkwd.imagAccumalate)/ngrid/nDataPerGridPt_) +
						" on average per data point");
		}
};

template<typename VTI, typename VTR>
void
FFTInterface::fft_data(
		VTI const & data,
		VTR & dataResult,
		int exponentSign)
{
	std::vector<int> gridDims = exponentSign > 0 ? gridDimsBKWD_ : gridDimsFWD_;
	int ngrid = 1;
	for ( auto di : gridDims )
		ngrid *= di;
	assert(data.size() == nDataPerGridPt_*ngrid);

	int threadID = omp_get_thread_num();

	this->allocate_this_thread(gridDims, nDataPerGridPt_);

	detail::ComplexConversion< std::complex<double>, typename VTI::value_type > converterFwd;

	auto fillBuffer = [&] () {
		auto ptr = reinterpret_cast<std::complex<double> * >(FFTBuffer_[threadID]);
		for (int ig = 0 ; ig < ngrid; ++ig)
			for (int id = 0 ; id < nDataPerGridPt_; ++id)
			{
				ptr[ id*ngrid + ig ] = converterFwd.convert( data [ ig * nDataPerGridPt_ + id ] );
				assert( ptr[ id*ngrid + ig ] == ptr[ id*ngrid + ig ] );
			}
	};

	if ( exponentSign > 0 )
	{
		fillBuffer();
		fftw_execute_dft(fftw3PlanBkwdKtoR_, FFTBuffer_[threadID], FFTBuffer_[threadID]);
		this->fill_result(ngrid, dataResult);
	}
	else
	{
		fillBuffer();
		fftw_execute_dft(fftw3PlanFowdRtoK_, FFTBuffer_[threadID], FFTBuffer_[threadID]);
		this->fill_result(ngrid, dataResult);
	}
}

template<typename VT>
void
FFTInterface::fft_interpolate(
		std::vector<int> const & gridDimsIn,
		std::vector<double> const & gridShiftIn,
		VT const & data,
		std::vector<int> const & gridDimsOut,
		std::vector<double> const & gridShiftOut,
		VT & dataResult,
		int nDataPerGridPt)
{
	typedef typename VT::value_type T;
	typedef typename detail::MakeComplex<typename VT::value_type>::type CT;

	int dim = gridDimsIn.size();
	assert(dim == gridDimsOut.size());
	int ngridIn = 1;
	int ngridOut = 1;
	for ( int id = 0 ; id < dim; ++id)
	{
		ngridIn *= gridDimsIn[id];
		ngridOut *= gridDimsOut[id];
	}
	assert( (ngridIn != 0) && (ngridOut != 0));

	// here we assume x major grid order
	std::vector<CT> intermedOut;
	this->plan_fft(gridDimsIn, nDataPerGridPt, -1, false, 1);
	this->fft_data(data, intermedOut, -1);

	for ( auto &d : intermedOut )
		d /= T(ngridIn);

	auto check_in_grid = [&] (std::vector<int> const & xyz) {
		for ( int id = 0 ; id < dim; ++id)
			if ( (xyz[id] <= -gridDimsOut[id]/2-gridDimsOut[id]%2) or (xyz[id] > gridDimsOut[id]/2) )
				return false;
		return true;
	};

	std::vector<CT> intermedIn(ngridOut*nDataPerGridPt, CT(0.0));
	std::vector<int> xyz(dim);
	for (int ic = 0 ; ic < ngridIn; ++ic )
	{
		FFTInterface::cnsq_to_xyz(ic, xyz, gridDimsIn, false);
		// go in to the frequency representation with + and - frequencies in the input mesh
		Algorithms::FFTInterface::inplace_to_freq(xyz, gridDimsIn);
		CT phase1 = CT(1);
		if (std::abs(*std::max_element(gridShiftIn	.begin(), gridShiftIn.end())) > 1e-6 )
		{
			auto sum = 0.0;
			for ( int id = 0; id < dim ; ++id)
				sum += (2*M_PI/gridDimsIn[id])*xyz[id]*gridShiftIn[id];
			phase1 = std::exp(CT(0,-sum));
		}

		// frequencies outside the output grid are dropped
		if ( not check_in_grid(xyz) )
			 continue;

		CT phase2 = CT(1);
		if (std::abs(*std::max_element(gridShiftOut.begin(), gridShiftOut.end())) > 1e-6 )
		{
			auto sum = 0.0;
			for ( int id = 0; id < dim ; ++id)
				sum += (2*M_PI/gridDimsOut[id])*xyz[id]*gridShiftOut[id];
			phase2 = std::exp(CT(0,sum));
		}

		// go in to the fftw3 representation in the positive frequency notation in the output mesh
		Algorithms::FFTInterface::freq_to_inplace(xyz, gridDimsOut);
		int cnsqNew = 0;
		Algorithms::FFTInterface::xyz_to_cnsq(cnsqNew, xyz, gridDimsOut, false);

		for ( int ib = 0 ; ib < nDataPerGridPt; ++ib )
			   intermedIn[cnsqNew*nDataPerGridPt+ib] =
					   intermedOut[ic*nDataPerGridPt+ib]*phase1*phase2;
	}

	// handle the Nyquist frequencies. The input data to the next FFT must be
	// balanced in frequencies. If the result was zero-padded, we simply add
	// the N/2 term to the -N/2 position and devide by 2, otherwise, we drop the
	// unbalanced term.
	for (int ic = 0 ; ic < ngridIn; ++ic )
	{
		FFTInterface::cnsq_to_xyz(ic, xyz, gridDimsIn, false);
		Algorithms::FFTInterface::inplace_to_freq(xyz, gridDimsIn);

		// Check which dimensions are on the non-symmetric edge while
		// there is zero padding to be used for leveling
		std::vector<int> copyDim;
		for ( int id = 0 ; id < dim; ++id)
			if ( ((gridDimsIn[id] % 2) == 0) and (gridDimsOut[id] > gridDimsIn[id]) )
				if ( xyz[id] == gridDimsIn[id]/2 )
					copyDim.push_back(id);

		if ( copyDim.empty() )
			continue;

		// Build the surface to be copied. While the below seems a bit opaque,
		// it solves problems in D >= 2 where the corner must be averaged between
		// 2^D points
		std::vector<std::vector<int>> indicesToCopy(1, xyz);
		int npts = std::pow(2,copyDim.size());
		indicesToCopy.reserve(npts);
		for ( int idc = 0 ; idc < copyDim.size(); ++idc )
		{
			std::vector<std::vector<int>> indicesDimMinusOne = indicesToCopy;
			for ( auto &d : indicesDimMinusOne )
				d[copyDim[idc]] = -gridDimsIn[copyDim[idc]]/2;
			indicesToCopy.insert(indicesToCopy.end(), indicesDimMinusOne.begin(), indicesDimMinusOne.end());
		}

		// Do the copying of Nyquist frequencies
		for ( int ip = 0 ; ip < npts; ++ip)
		{
			// recompute phases for the copied Nyquist frequency
			CT phase1 = CT(1);
			if (std::abs(*std::max_element(gridShiftIn.begin(), gridShiftIn.end())) > 1e-6 )
			{
				auto sum = 0.0;
				for ( int id = 0; id < dim ; ++id)
					sum += (2*M_PI/gridDimsIn[id])*indicesToCopy[ip][id]*gridShiftIn[id];
				phase1 = std::exp(CT(0,-sum));
			}
			CT phase2 = CT(1);
			if (std::abs(*std::max_element(gridShiftOut.begin(), gridShiftOut.end())) > 1e-6 )
			{
				auto sum = 0.0;
				for ( int id = 0; id < dim ; ++id)
					sum += (2*M_PI/gridDimsOut[id])*indicesToCopy[ip][id]*gridShiftOut[id];
				phase2 = std::exp(CT(0,sum));
			}

			int cnsq;
			FFTInterface::freq_to_inplace(indicesToCopy[ip], gridDimsOut);
			FFTInterface::xyz_to_cnsq(cnsq, indicesToCopy[ip], gridDimsOut, false);
			for ( int ib = 0 ; ib < nDataPerGridPt; ++ib )
				intermedIn[cnsq*nDataPerGridPt+ib] =
						   intermedOut[ic*nDataPerGridPt+ib]*phase1*phase2 / CT(npts);
		}

		// We have to balance the output by dropping the Nyquist frequency if there is no space
		for ( int id = 0 ; id < dim; ++id)
			if ( ((gridDimsOut[id] % 2) == 0) and (gridDimsOut[id] < gridDimsIn[id]) )
			{
				FFTInterface::cnsq_to_xyz(ic, xyz, gridDimsOut, false);
				if ( xyz[id] == gridDimsOut[id]/2 )
					for ( int ib = 0 ; ib < nDataPerGridPt; ++ib )
						intermedIn[ic*nDataPerGridPt+ib] = 0.0;
			}
	}
	this->plan_fft(gridDimsOut, nDataPerGridPt, 1, false, 1);
	this->fft_data(intermedIn, dataResult, 1);
}

template<typename T>
void
FFTInterface::fft_hessian(
		std::vector<double> const & latticeMatrix,
		std::vector<int> const & grid,
		std::vector< T > const & data,
		std::vector< T > & hessianOfData,
		int nDataPerGridPt)
{
	typedef typename detail::MakeComplex<T>::type CT;

	int dim = grid.size();
	int nG = 1;
	for ( int id = 0 ; id < dim; ++id)
		nG *= grid[id];
	assert(nG != 0);
	assert(latticeMatrix.size() == dim*dim);
	assert(data.size() == nG*nDataPerGridPt);

	// here we assume x major grid order
	std::vector<CT> intermedOut;
	this->plan_fft(grid, nDataPerGridPt, -1, false, 1);
	this->fft_data(data, intermedOut, -1);

	std::vector<CT> fftHessian(nDataPerGridPt*nG*dim*dim, CT(0));
	std::vector<int> xyz(dim);
	std::vector<double> rvec(dim);
	for (int ig = 0 ; ig < nG; ++ig )
	{
		std::fill(rvec.begin(), rvec.end(), 0.0);
		cnsq_to_xyz(ig, xyz, grid, false);
		inplace_to_freq(xyz, grid);
		for ( int i = 0 ; i < dim ; ++i)
			for ( int j = 0 ; j < dim ; ++j)
				rvec[i] += latticeMatrix[i*dim + j]*xyz[j];

		for ( int id = 0 ; id < nDataPerGridPt ; ++id)
			for ( int i = 0 ; i < dim ; ++i)
				for ( int j = 0 ; j < dim ; ++j)
				{
					int csq = j + dim*(i + dim*(id + nDataPerGridPt*ig));
					fftHessian[csq] = -intermedOut[ig*nDataPerGridPt+id]*rvec[i]*rvec[j]/CT(nG);
				}
	}
	this->plan_fft(grid, nDataPerGridPt*dim*dim, 1, false, 1);
	this->fft_data(fftHessian, hessianOfData, 1);
}

} /* namespace Algorithms */
} /* namespace elephon */
