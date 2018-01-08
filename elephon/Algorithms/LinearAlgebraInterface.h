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
#ifdef USE_MKL

#define MKL_INT int
#define MKL_Complex16 std::complex<double>
#define MKL_Complex8 std::complex<float>
#include <mkl_types.h>
#include <mkl_lapacke.h>
#include <mkl.h>

#define CBLAS_CMPLDBL_PTR void
#define CBLAS_CMPLFLT_PTR void

#elif USE_ATLAS
extern "C"
{
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include <cblas.h>
#include <lapacke.h>

#define CBLAS_CMPLDBL_PTR void
#define CBLAS_CMPLFLT_PTR void
}
#else
extern "C"
{
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#define CBLAS_CMPLDBL_PTR double
#define CBLAS_CMPLFLT_PTR float
#include <cblas.h>
#include <lapacke.h>
}
#endif
#include "Auxillary/AlignedVector.h"
#include <vector>

namespace elephon
{
namespace Algorithms
{

namespace detail
{
template<typename T>
struct ComplexTypeTrait
{
	typedef T type;
};

template<typename T>
struct ComplexTypeTrait< std::complex<T> >
{
	typedef T type;
};
}

/**
 * This linear algebra interface is a total mess. Passes tests and all but is really not in a good shape ...
 * TODO replace it with some usable library.
 *	Candidates for a replacement: FLENS, The Armadillo C++ linear algebra library
 *
 */
class LinearAlgebraInterface
{
public:
	//Here come the high level routines
	template<typename T>
	void pseudo_inverse(std::vector<T> A, int n, int m,
			std::vector<T> & pinvA, double cutoff = 1e-10);

	template<typename T>
	void null_space(std::vector<T> A, int n, int m,
			int & kerDim, std::vector<T> & nullA, double tol = 1e-5);

	template<class C>
	void inverse(C A, C & invA );

	template<typename T>
	void diagonalize_hermitian(bool data_upper, bool comEV,
			std::vector< std::complex<T> > matrix,
			std::vector< std::complex<T> > & eigenvectors,
			std::vector<T> & eigenvalues );

	template<typename VT>
	void matrix_matrix_prod(VT const & A,
			VT const & B,
			VT & ATimesB, int m, int n) const;

	template<typename T>
	void svd(std::vector<T> A, int m, int n,
			std::vector<T> & U,
			std::vector<T> & VT,
			std::vector< typename detail::ComplexTypeTrait< std::complex<T> >::type > & sv);

	template<typename VT>
	void transpose_square_matrix(VT & A) const;

	template<typename VT>
	void conjugate_square_matrix(VT & A) const;


	//Here come the low level routines

	/**
	 * GEMV  performs one of the matrix-vector operations
	 *
	 * y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
	 *
	 *	where alpha and beta are scalars, x and y are vectors and A is an
	 *	m by n matrix.
	 */
	void call_gemv(
			char TRANS,
			int M,
			int N,
			std::complex<float> ALPHA,
			std::complex<float> const *	A,
			int LDA,
			std::complex<float> const * X,
			int INCX,
			std::complex<float> BETA,
			std::complex<float> * Y,
			int  INCY) const;

	std::complex<float> call_dotu(
			int n,
			std::complex<float> * dx, int incx,
			std::complex<float> * dy, int incy) const;

	std::complex<double> call_dotu(
			int n,
			std::complex<double> * dx, int incx,
			std::complex<double> * dy, int incy) const;

	int call_gemm(
			char transA, char transB,
	        int m, int n, int k, double alpha, double const * A, int lda,
	        double const * B, int ldb, double beta,double * C,int ldc) const;

	int call_gemm(
			char transA, char transB,
	        int m, int n, int k, std::complex<float> alpha, std::complex<float> const * A, int lda,
			std::complex<float> const * B, int ldb, std::complex<float> beta,std::complex<float> * C,int ldc) const;

	int call_gemm(
			char transA, char transB,
	        int m, int n, int k, std::complex<double> alpha, std::complex<double> const * A, int lda,
			std::complex<double> const * B, int ldb, std::complex<double> beta,std::complex<double> * C,int ldc) const;

	int call_getri( int matrix_order, int n, double * a, int lda,
			const int * ipiv, double * work, int lwork);

	int call_getrf( int matrix_order, int m, int n, double * a, int lda, int * ipiv );

	int call_heev( int matrix_order, char jobz, char uplo,
			   int n, std::complex<double> * a,
			   int lda, double * w,
			   std::complex<double> * work, int lwork,
			   double * rwork) const;

	int call_heevr(  int matrix_order, char jobz, char uplo,
			char range,  int n, std::complex<double> * a, int lda, double vl, double vu,
			int il, int iu, double abstol, int * m, double * w, std::complex<double> * z, int ldz, int * isuppz,
			 std::complex<double> * work, int lwork,
			double * rwork, int lrwork, int * iwork, int liwork) const;

	int call_gesdd(
			int matrix_order, char jobz, int m,
			int n, double* a, int lda,
	        double* s, double* u, int ldu,
	        double* vt, int ldvt, double* work,
			double * rwork,
			int lwork, int* iwork  );

	int call_syev(
			char jobz, char uplo, int n, double * a, int lda, double * w,
			int matrix_layout = CblasRowMajor);

	int call_syev(
			char jobz, char uplo, int n, float * a, int lda, float * w,
			int matrix_layout = CblasRowMajor);

	void clear_buffer();
private:

	std::vector<int> IPIV_;

	Auxillary::alignedvector::CV workbuffer_;

	Auxillary::alignedvector::CV rWork_;

	template<class C>
	int square_matrix_dim(C const & A) const;

	void check_library_info(int returnCode, std::string const & calledByWhat ) const;

	CBLAS_TRANSPOSE ctoen(char const & c) const;
};

} /* namespace Algorithms */
} /* namespace elephon */

#include "Algorithms/LinearAlgebraInterface.hpp"
#endif /* ELEPHON_ALGORITHMS_LINEARALGEBRAINTERFACE_H_ */
