/*	This file FFTInterface.h is part of elephon.
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

#ifndef ELEPHON_ALGORITHMS_FFTINTERFACE_H_
#define ELEPHON_ALGORITHMS_FFTINTERFACE_H_

#include "LatticeStructure/LatticeModule.h"
#include <vector>
#include <complex>
#include <fftw3.h>

namespace elephon
{
namespace Algorithms
{

/**
 * Interface to the fftw3 library that implements the fast Fourier transform (FFT).
 *
 * This interface not only provides the high level methods to compute the FFT is also buffers
 * the storage used by the FFTW library for fast calculation of multiple similar FFTs.
 * This class is thread-safe for OMP threaded applications.
 * NOTE: For other methods, this interface is NOT thread safe.
 */
class FFTInterface
{
public:

	/**
	 * Clears the internal storage.
	 */
	~FFTInterface();

	/**
	 * Compute a regular grid Fourier transform from a packed storage scheme.
	 *
	 * @tparam VI	must be vectors with possibly custom allocators. Linear consecutive memory layout of the data
	 * 				is necessary as required by the underlying fftw library..
	 * @tparam VR	must be vectors with possibly custom allocators. Linear consecutive memory layout of the data
	 * 				is necessary as required by the underlying fftw library..
	 *
	 * @param[in] mapFFTCoeff			A list of 3N elements gx1,gy1,gz1,gx2,...,gzN specifying the G vector of
	 * 								the respective data values in \ref sparseInputData
	 * @param[in] gridDimsInputData		The 3 element vector (NGX, NGY, NGZ) with the regular grid dimensions that the \p sparseInputData refers to.
	 * @param[in] sparseInputData		The data values V1, V2 ... VN corresponding to the G vectors in \p mapFFTCoeff
	 * @param[in] exponentSign			The sign of the complex exponent.
	 * @param[out] dataResult			A list of transformed data on the grid and data layout as specified in the call to plan_fft()
	 */
	template<typename VI, typename VR>
	void fft_sparse_data(
			std::vector<int> const & mapFFTCoeff,
			std::vector<int> const & gridDimsInputData,
			VI const & sparseInputData,
			int exponentSign,
			VR & dataResult);

	void plan_fft(
			std::vector<int> const & gridDimsData,
			int nDataPerGridPt = 1,
			int exponentSign = 0,
			bool dataLayoutRowMajor = false,
			int hintHowOften = 1);

	/**
	 * Perform the Fourier transform of data.
	 *
	 * @param data				the data in the layout provided when plan_fft was called. For each
	 * 							grid point, nDataPerGridPt are assumed to be in consecutive places.
	 * 							If the grid is e.g. (x,y,z) the data layout is for column major
	 * 								consq index = idata + nDataPerGridPt*(ix + x*(iy + y*(iz) ) )
  	 *							nDataPerGridPt is set by the planner with a call to plan_fft
	 * @param dataResult		vector that will be resized to fit the Fourier transformed data.
	 * @param exponentSign		if the transformation is forward or backwards. Must not be zero.
	 */
	template<typename VTI, typename VTR>
	void fft_data(
			VTI const & data,
			VTR & dataResult,
			int exponentSign);

	/**
	 * Perform a band-limited Fourier interpolation of data and place it in dataResult.
	 *
	 * @param gridDimsIn		Input grid with integer indicating the size, e.g. {5 5}
	 * @param gridShiftIn		Grid shift in units of the grid devision. Must have the same size as \p gridDimsIn
	 * @param data				The input data is assumed on the grid \p gridDimsIn shifted by \p gridShiftIn
	 * 							in dimension nDataPerGridPt, x[, y[,...]] laid out in Column major or Fortran order.
	 * @param gridDimsOut		The output grid. Must be of the same size as \p gridDimsIn
	 * @param gridShiftOut		Grid shift in units of the grid devision. Must have the same size as \p gridDimsIn
	 * @param dataResult		The interpolated \p data. Same data order as input. Will be resized to fit the output.
	 * @param nDataPerGridPt	The number of data elements per grid point.
	 */
	template<typename VT>
	void fft_interpolate(
			std::vector<int> const & gridDimsIn,
			std::vector<double> const & gridShiftIn,
			VT const & data,
			std::vector<int> const & gridDimsOut,
			std::vector<double> const & gridShiftOut,
			VT & dataResult,
			int nDataPerGridPt);

	/**
	 * Compute the second derivative matrix of a data field.
	 *
	 * @param latticeMatrix		The basis of the reciprocal space where the field is defined.
	 * 							Must be a vector of size dimension^2, where 'dimension' is the
	 * 							size of the vector of parameter \p grid with the basis vectors of size dimension
	 * 							in fastest running order. Thus in 3D the first vector must be in position 0, 1 and 2.
	 * @param grid				The number of elements in each dimension x[,y[,z[...]].
	 * @param data				Vector of size (Prod grid) times \p nDataPerGridPt. Input data, assumed in Column major grid
	 * 							order such that each block of \p nDataPerGridPt data points are layed out in x fastest order.
	 * @param hessianOfData		Output data. On output will be a vector of size (Prod grid) times \p nDataPerGridPt times
	 * 							dim*dim. Will be resized to hold the data. The layout is
	 * 							[xi,xj,d,grid] where xi and xj enumare the dimension, d the data points per grid value and
	 * 							grid is the the grid in Column major order.
	 * @param nDataPerGridPt	Number of data elements per grid point.
	 */
	template<typename T>
	void fft_hessian(
			std::vector<double> const& latticeMatrix,
			std::vector<int> const & grid,
			std::vector< T > const & data,
			std::vector< T > & hessianOfData,
			int nDataPerGridPt);

	/**
	 * Release internal stored objects, such that a new plan_fft is necessary.
	 */
	void clear_storadge();


	static void inplace_to_freq(int &x, int &y, int &z, int dx, int dy, int dz);

	static void inplace_to_freq(std::vector<int> & g, std::vector<int> const & d);

	static void freq_to_inplace(int &x, int &y, int &z, int dx, int dy, int dz);

	static void freq_to_inplace(std::vector<int> & g, std::vector<int> const & d);

	static void cnsq_to_xyz(int cnsq,
			std::vector<int> & xyz,
			std::vector<int> const& grid,
			bool dataLayoutRowMajor = false);

	static void xyz_to_cnsq(int &cnsq,
			std::vector<int> const& xyz,
			std::vector<int> const& grid,
			bool dataLayoutRowMajor = false);
private:

	//do not implement
	FFTInterface(FFTInterface const & );
	FFTInterface * operator= (FFTInterface);

	int numThreadsMax_ = 1;

	std::vector<int> nBuff_;

	std::vector<int> gridDimsFWD_;

	std::vector<int> gridDimsBKWD_;

	bool dataLayoutRowMajor_ = false;

	int nDataPerGridPt_ = 1;

	fftw_complex ** FFTBuffer_ = nullptr;

	fftw_plan_s * fftw3PlanFowdRtoK_ = nullptr;

	fftw_plan_s * fftw3PlanBkwdKtoR_ = nullptr;

	void allocate_internal();

	void allocate_this_thread(
			std::vector<int> const & gridDims,
			int nDataPerGridPt);

	template< typename VTR>
	void fill_result(int ngrid, VTR & dataResult );

	void plan_fft_local(
			std::vector<int> const & gridDims,
			int exponentSign,
			int nDataPerGridPt,
			int hintHowOften,
			fftw_plan_s * & toAlloc,
			bool dataLayoutRowMajor );
};

} /* namespace Algorithms */
} /* namespace elephon */

#include "Algorithms/FFTInterface.hpp"
#endif /* ELEPHON_ALGORITHMS_FFTINTERFACE_H_ */
