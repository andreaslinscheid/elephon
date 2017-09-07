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

class FFTInterface
{
public:
	FFTInterface();

	~FFTInterface();

	template<typename TI, typename TR>
	void fft_sparse_data(
			std::vector<int> const & mapFFTCoeff,
			std::vector<int> const & gridDimsInputData,
			std::vector< TI > const & sparseInputData,
			int nDataPerGridPt,
			int exponentSign,
			std::vector< TR > & dataResult,
			std::vector<int> const & gridDimsOutputData,
			bool dataLayoutRowMajor = false,
			int hintHowOften = 1);

	template<typename TI, typename TR>
	void fft_data(
			std::vector<int> const & gridDimsData,
			std::vector< TI > const & data,
			std::vector< TR > & dataResult,
			int nDataPerGridPt,
			int exponentSign,
			bool dataLayoutRowMajor = false,
			int hintHowOften = 1);

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
	template<typename T>
	void fft_interpolate(
			std::vector<int> const & gridDimsIn,
			std::vector<double> const & gridShiftIn,
			std::vector< T > const & data,
			std::vector<int> const & gridDimsOut,
			std::vector<double> const & gridShiftOut,
			std::vector< T > & dataResult,
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

	int nBuff_ = 0;

	std::vector<int> gridDimsFWD_;

	std::vector<int> gridDimsBKWD_;

	fftw_complex * FFTBuffer_ = nullptr;

	fftw_plan_s * fftw3PlanFowdRtoK_ = nullptr;

	fftw_plan_s * fftw3PlanBkwdKtoR_ = nullptr;

	template< typename TR>
	void allocate(
			std::vector<int> const & gridDims,
			int nDataPerGridPt,
			std::vector< TR > & dataResult,
			bool & re_plan);

	template< typename TR>
	void fill_result(int ngrid, int nDataPerGridPt, std::vector<TR> & dataResult );

	void plan_fft(
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
