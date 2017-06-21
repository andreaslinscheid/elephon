/*	This file LinearAlgebraInterface.h is part of elephon.
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
 *  Created on: Jun 7, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_ALGORITHMS_LINEARALGEBRAINTERFACE_H_
#define ELEPHON_ALGORITHMS_LINEARALGEBRAINTERFACE_H_

#include <complex>
#include <vector>

namespace elephon
{
namespace Algorithms
{

class LinearAlgebraInterface
{
public:
	//Here come the high level routines
	template<typename T>
	void pseudo_inverse(std::vector<T> A, int n, int m,
			std::vector<T> & pinvA, double cutoff = 1e-10);

	template<typename T>
	void inverse(std::vector<T> A, std::vector<T> & invA );

	template<typename T>
	void diagonalize_hermitian(bool data_upper, bool comEV,
			std::vector< std::complex<T> > matrix,
			std::vector< std::complex<T> > & eigenvectors,
			std::vector<T> & eigenvalues );

	template<typename T>
	void matrix_matrix_prod(std::vector<T> const & A,
			std::vector<T> const & B,
			std::vector<T> & ATimesB, int m, int n) const;

	//Here come the low level routines

	int call_gemm(
			char transA, char transB,
	        int m, int n, int k, double alpha, double const * A, int lda,
	        double const * B, int ldb, double beta,double * C,int ldc) const;

	int call_getri( int matrix_order, int n, double * a, int lda,
			const int * ipiv, double * work, int lwork);

	int call_getrf( int matrix_order, int m, int n, double * a, int lda, int * ipiv );

	int call_heev( int matrix_order, char jobz, char uplo,
			   int n, std::complex<double> * a,
			   int lda, double * w,
			   std::complex<double> * work, int lwork,
			   double * rwork) const;

	int call_gesdd(
			int matrix_order, char jobz, int m,
			int n, double* a, int lda,
	        double* s, double* u, int ldu,
	        double* vt, int ldvt, double* work,
			double * rwork,
			int lwork, int* iwork  );

	void clear_buffer();
private:

	std::vector<int> IPIV_;

	std::vector<char> workbuffer_;

	std::vector<char> rWork_;

	template<class C>
	int square_matrix_dim(C const & A) const;

	void check_library_info(int returnCode, std::string const & calledByWhat ) const;
};

} /* namespace Algorithms */
} /* namespace elephon */

#include "Algorithms/LinearAlgebraInterface.hpp"
#endif /* ELEPHON_ALGORITHMS_LINEARALGEBRAINTERFACE_H_ */
