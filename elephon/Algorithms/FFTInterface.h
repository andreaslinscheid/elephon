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

	template<typename T>
	void fft_interpolate(
	                std::vector<int> const & gridDimsIn,
	                std::vector< T > const & data,
	                std::vector<int> const & gridDimsOut,
	                std::vector< T > & dataResult,
	                int nDataPerGridPt);


	static void inplace_to_freq(int &x, int &y, int &z, int dx, int dy, int dz);

	static void inplace_to_freq(std::vector<int> & g, std::vector<int> const & d);

	static void freq_to_inplace(int &x, int &y, int &z, int dx, int dy, int dz);

	static void freq_to_inplace(std::vector<int> & g, std::vector<int> const & d);
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
