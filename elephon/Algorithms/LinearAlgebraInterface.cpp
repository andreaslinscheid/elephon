/*	This file LinearAlgebraInterface.cpp is part of elephon.
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

#include "LinearAlgebraInterface.h"
#include <cmath>
#include <assert.h>
#include <stdexcept>

namespace elephon
{
namespace Algorithms
{

void
LinearAlgebraInterface::clear_buffer()
{
	IPIV_.clear();
	workbuffer_.clear();
	rWork_.clear();
}

void
LinearAlgebraInterface::call_gemv(
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
		int INCY) const
{
	auto cblasTrans = CblasNoTrans;
	if ( (TRANS == 'T') or (TRANS == 't') )
		cblasTrans = CblasTrans;
	if ( (TRANS == 'C') or (TRANS == 'c') )
		cblasTrans = CblasConjTrans;
	cblas_cgemv(CblasRowMajor, cblasTrans, M, N,
			reinterpret_cast<CBLAS_CMPLFLT_PTR*>(&ALPHA),
			reinterpret_cast<const CBLAS_CMPLFLT_PTR*>(A), LDA,
			reinterpret_cast<const CBLAS_CMPLFLT_PTR*>(X), INCX,
			reinterpret_cast<CBLAS_CMPLFLT_PTR*>(&BETA), reinterpret_cast<CBLAS_CMPLFLT_PTR *>(Y), INCY);
}

std::complex<float>
LinearAlgebraInterface::call_dotu(
		int n,
		std::complex<float> * dx, int incx,
		std::complex<float> * dy, int incy) const
{
	std::complex<float> result;
#if defined USE_MKL || defined USE_ATLAS
	cblas_cdotu_sub(n, reinterpret_cast<const void*>(dx), incx,
					   reinterpret_cast<const void*>(dy), incy,
					   reinterpret_cast<void*>(&result));
#else // openblas
	cblas_cdotu_sub(n, reinterpret_cast<const CBLAS_CMPLFLT_PTR*>(dx), incx,
					   reinterpret_cast<const CBLAS_CMPLFLT_PTR*>(dy), incy,
					   reinterpret_cast<openblas_complex_float*>(&result));
#endif
	return result;
}

std::complex<double>
LinearAlgebraInterface::call_dotu(
		int n,
		std::complex<double> * dx, int incx,
		std::complex<double> * dy, int incy) const
{
	std::complex<double> result;
#if defined USE_MKL || defined USE_ATLAS
	cblas_zdotu_sub(n, reinterpret_cast<const void*>(dx), incx,
					   reinterpret_cast<const void*>(dy), incy,
					   reinterpret_cast<void*>(&result));
#else// openblas
	cblas_zdotu_sub(n, reinterpret_cast<const CBLAS_CMPLDBL_PTR*>(dx), incx,
					   reinterpret_cast<const CBLAS_CMPLDBL_PTR*>(dy), incy,
					   reinterpret_cast<openblas_complex_double*>(&result));
#endif
	return result;
}

int
LinearAlgebraInterface::call_gemm(
		char transA, char transB,
        int m, int n, int k, double alpha, double const * A, int lda,
        double const * B, int ldb, double beta,double * C,int ldc) const
{
	cblas_dgemm(CblasRowMajor, ctoen(transA), ctoen(transB),
			m, n, k,
			alpha, A, lda,
			B, ldb,
			beta, C, ldc);
	return 0;
}

CBLAS_TRANSPOSE
LinearAlgebraInterface::ctoen(char const & c) const
{
	switch (c) {
		case 'c':
		case 'C':
			return CblasConjTrans;
		case 'n':
		case 'N':
			return CblasNoTrans;
		case 't':
		case 'T':
			return CblasTrans;
		default:
			throw std::runtime_error("Invalid char passed to gemm");
			break;
	}
	return CblasNoTrans;
}

int
LinearAlgebraInterface::call_gemm(
		char transA, char transB,
        int m, int n, int k, std::complex<float> alpha, std::complex<float> const * A, int lda,
		std::complex<float> const * B, int ldb, std::complex<float> beta,std::complex<float> * C,int ldc) const
{
	cblas_cgemm(CblasRowMajor, ctoen(transA), ctoen(transB),
			m, n, k,
			reinterpret_cast<CBLAS_CMPLFLT_PTR*>(&alpha), reinterpret_cast<const CBLAS_CMPLFLT_PTR*>(A), lda,
			reinterpret_cast<const CBLAS_CMPLFLT_PTR*>(B), ldb,
			reinterpret_cast<CBLAS_CMPLFLT_PTR*>(&beta), reinterpret_cast<CBLAS_CMPLFLT_PTR*>(C), ldc);
	return 0;
}

int
LinearAlgebraInterface::call_gemm(
		char transA, char transB,
        int m, int n, int k, std::complex<double> alpha, std::complex<double> const * A, int lda,
		std::complex<double> const * B, int ldb, std::complex<double> beta,std::complex<double> * C,int ldc) const
{
	cblas_zgemm(CblasRowMajor, ctoen(transA), ctoen(transB),
			m, n, k,
			reinterpret_cast<CBLAS_CMPLDBL_PTR*>(&alpha), reinterpret_cast<const CBLAS_CMPLDBL_PTR*>(A), lda,
			reinterpret_cast<const CBLAS_CMPLDBL_PTR*>(B), ldb,
			reinterpret_cast<CBLAS_CMPLDBL_PTR*>(&beta), reinterpret_cast<CBLAS_CMPLDBL_PTR*>(C), ldc);
	return 0;
}

int
LinearAlgebraInterface::call_getri(
		int matrix_order, int n, double * a, int lda,
		const int * ipiv, double * work, int lwork)
{
	return LAPACKE_dgetri_work(matrix_order,n,a,lda,ipiv,work,lwork);
}

int
LinearAlgebraInterface::call_getrf( int matrix_order, int m, int n, double * a, int lda, int * ipiv )
{
	return LAPACKE_dgetrf_work(matrix_order,m,n,a,lda,ipiv);
}

int
LinearAlgebraInterface::call_heev( int matrix_order, char jobz, char uplo,
		   int n, std::complex<double> * a,
		   int lda, double * w,
		   std::complex<double> * work, int lwork,
		   double * rwork) const
{
	return LAPACKE_zheev_work(matrix_order,jobz, uplo, n, a, lda,w,work,lwork,rwork);
}

int
LinearAlgebraInterface::call_heevr( int matrix_order, char jobz, char uplo,
		char range,  int n, std::complex<double> * a, int lda, double vl, double vu,
		int il, int iu, double abstol, int * m, double * w, std::complex<double> * z, int ldz, int * isuppz,
		 std::complex<double> * work, int lwork,
		double * rwork, int lrwork, int * iwork, int liwork) const
{
	return LAPACKE_zheevr_work(matrix_order, jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z,
			ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork);
}

int
LinearAlgebraInterface::call_gesdd(
		int matrix_order, char jobz, int m,
		int n, double* a, int lda,
        double* s, double* u, int ldu,
        double* vt, int ldvt, double* work,
		double * rwork, //not used in the real version - only in complex
		int lwork, int* iwork  )
{
	return LAPACKE_dgesdd_work( matrix_order, jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork );
}

int
LinearAlgebraInterface::call_syev(char jobz, char uplo, int n, double * a, int lda, double * w, int matrix_layout)
{
	return LAPACKE_dsyev(matrix_layout, jobz, uplo, n, a, lda, w);
}

int
LinearAlgebraInterface::call_syev(char jobz, char uplo, int n, float * a, int lda, float * w, int matrix_layout)
{
	return LAPACKE_ssyev(matrix_layout, jobz, uplo, n, a, lda, w);
}

void
LinearAlgebraInterface::check_library_info(int returnCode, std::string const & calledByWhat ) const
{
	if ( returnCode == 0 )
		return;

	if ( returnCode < 0 )
		throw std::runtime_error( std::string("input ")+calledByWhat+" argument # "
				+std::to_string(-returnCode)+"  had an illegal value !");

	if ( returnCode > 0 )
		throw std::runtime_error( std::string("Algorithm ")+calledByWhat+"failed to converge:"
				+std::to_string(-returnCode)+" !");
}

} /* namespace Algorithms */
} /* namespace elephon */
