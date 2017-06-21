/*	This file LinearAlgebraInterface.hpp is part of elephon.
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
 *  Created on: Jun 8, 2017
 *      Author: A. Linscheid
 */

#include "LinearAlgebraInterface.h"
#ifdef USE_MKL

#define MKL_INT int
#define MKL_Complex16 std::complex<double>
#define MKL_Complex8 std::complex<float>
#include <mkl_types.h>
#include <mkl_lapacke.h>
#include <mkl.h>

#else
extern "C"
{
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include <cblas.h>
#include <lapacke.h>
}
#endif
#include <cmath>
#include <assert.h>
#include <stdexcept>

#include <iostream>

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

template<typename T>
void
LinearAlgebraInterface::pseudo_inverse(std::vector<T> A, int m, int n,
		std::vector<T> & pinvA, double cutoff)
{
	//This routine needs to buffer a few matrices.
	//In a kinda hacky fashion, we put them in the workspace, too.
	typedef typename detail::ComplexTypeTrait<T>::type realT;

	//
	//	First compute the SVD
	//
	pinvA = std::move( A );
	int ldu =  m;
	int ldvt = n;
	int lda = n;
	int info;
	IPIV_.resize( std::max(1, 8 * std::min(m, n)) );
	T work_query;
	info = this->call_gesdd(LAPACK_ROW_MAJOR,'A',m,n,NULL,lda,NULL,NULL,ldu,NULL,ldvt,&work_query,NULL,-1,IPIV_.data());
	check_library_info(info,"gesdd");

	//TODO here we should worry about alignment!
	int lwork = static_cast<int>( std::real(work_query) );
	int nElem = lwork + (m*m)/*U matrix*/ + (n*n)/*V^T matrix*/ + (n*m) /*U * \Sigma+*/;
	int nElemChar = sizeof(T)*nElem + sizeof(realT)*std::min(m, n)/*Singular values*/;
	workbuffer_.resize( nElemChar );

	//We put work first, u after that and v^t and finally the singular values
	T * work = reinterpret_cast< T * >( &workbuffer_[0] );
	T * u = reinterpret_cast< T * >( &workbuffer_[ sizeof(T)*lwork] );
	T * vt = reinterpret_cast< T * >( &workbuffer_[ sizeof(T)*(lwork+m*m) ] );
	T * uSigmaP = reinterpret_cast< T * >( &workbuffer_[ sizeof(T)*(lwork+m*m+n*n)] );
	realT * s = reinterpret_cast< realT * >( &workbuffer_[ sizeof(T)*nElem ] );

	if ( sizeof(realT) != sizeof(T) )
		//We are running the complex version
		rWork_.resize( sizeof(realT)*(std::max(1,std::min(m,n)*std::max(5*std::min(m,n)+7,2*std::max(m,n)+2*std::min(m,n)+1))) );
		//The above formula (mess) is from https://software.intel.com/en-us/mkl-developer-reference-fortran-gesdd

	//The complex version needs a parameter rwork
	realT * rwork =  reinterpret_cast< realT * >( rWork_.data() );

	//Left singular vectors (U) are stored column wise, right singular vectors (vT) are stored row wise
	this->call_gesdd( LAPACK_ROW_MAJOR, 'A', m, n, pinvA.data(), lda, s, u, ldu, vt, ldvt, work, rwork, lwork, IPIV_.data()  );

	//
	//	Compute (v^T)^T * \Sigma+ * U^T = v * ( U * \Sigma+^T )^T
	//
	realT relCut = cutoff*s[0];//the singular values must be returned sorted by lapack with the first one the largest

	//Make \Sigma into \Sigma+
	for ( int i = 0; i < std::min(n,m); ++i )
		s[i] = s[i] < relCut ? 0 : 1.0/s[i];

	//Multiplying first ( U * \Sigma+^T )
	std::fill(uSigmaP,uSigmaP+n*m,0.0);
	for ( int i = 0 ; i < m; ++i)
		for ( int j = 0 ;j < std::min(n,m); ++j)
			uSigmaP[i*n+j] += u[i*m+j]*s[j];

	//Now compute A+ = (v^T)^T * ( U * \Sigma+^T )^T
	this->call_gemm('c','c',n,m,n,T(1.0),vt,n,uSigmaP,n,T(0.0),pinvA.data(),m);
}

template<typename T>
void
LinearAlgebraInterface::matrix_matrix_prod(std::vector<T> const & A,
		std::vector<T> const & B,
		std::vector<T> & ATimesB, int m, int n) const
{
	int k = int(A.size())/m;
	assert( k == int(B.size())/n );
	ATimesB.resize( m * n );
	this->call_gemm('n','n',m,n,k,T(1.0),A.data(),k,B.data(),n,T(0.0),ATimesB.data(),n);
}

template<typename T>
void
LinearAlgebraInterface::inverse(std::vector<T> A, std::vector<T> & invA )
{
	invA = std::move(A);
	int dim = this->square_matrix_dim(invA);
	int info;
	if ( IPIV_.size() != static_cast<size_t>(dim) )
	{
		IPIV_.assign( dim, 0);
		double work_query;
		info = this->call_getri(LAPACK_ROW_MAJOR,dim,NULL,dim,NULL,&work_query,-1);
		check_library_info(info,"i 1 getri");
		workbuffer_.resize( sizeof(T)*static_cast<int>( work_query ) );
	}

	auto * dptr = invA.data();
	int * ipptr = IPIV_.data();
	auto * work = reinterpret_cast<T*>(workbuffer_.data());
	int lwork = static_cast<int>(workbuffer_.size())/sizeof(T);

	// L U factorization
	info = this->call_getrf(LAPACK_ROW_MAJOR,dim,dim,dptr,dim,ipptr);
	this->check_library_info(info,"i getrf");

	//compute the inverse
	info = this->call_getri(LAPACK_ROW_MAJOR,dim,dptr,dim,ipptr,work,lwork);
	check_library_info(info,"i 2 getri");
}

template<typename T>
void
LinearAlgebraInterface::diagonalize_hermitian(
		bool data_upper, bool comEV,
		std::vector< std::complex<T> > matrix,
		std::vector< std::complex<T> > & eigenvectors,
		std::vector<T> & eigenvalues )
{
	char v = ( comEV ?  'v' : 'n' );

	char dataLoc = ( data_upper ? 'U' : 'L' );

	int dim = this->square_matrix_dim( matrix );

	int rworksize = ( 3*dim-2 > 1 ? 3 * dim - 2 : 1 ) ;
	rWork_.resize(sizeof(T)*rworksize);

	T work_query;
	this->call_heev(
			LAPACK_ROW_MAJOR,
			v, dataLoc,
			dim, matrix, dim,
			eigenvalues.data(),&work_query,-1,
			reinterpret_cast< T* >(rWork_.data()));
	size_t optimal_size = static_cast<size_t>( work_query.real() ) ;
	workbuffer_.resize( sizeof(std::complex<T>)*optimal_size );

	int lwork = static_cast<int>(optimal_size);
	int info = this->call_heev(
			LAPACK_ROW_MAJOR,
			v, dataLoc,
			dim, matrix, dim,
			eigenvalues.data(),
			reinterpret_cast< std::complex<T>* >(workbuffer_.data()),
			lwork,
			reinterpret_cast< T* >(rWork_.data()));

	this->check_library_info(info, "heev diag");

	//Ensure phase convention on eigenvectors (real part of the last component is non-negative)
	for ( int ib = 0 ; ib < dim ; ++ib )
		if ( matrix[(dim-1)*dim + ib].real() < 0 )
			for ( int ibp = 0 ; ibp < dim ; ++ibp )
				matrix[ibp*dim + ib] *= std::complex<T>(-1.0);
}

template<class C>
int
LinearAlgebraInterface::square_matrix_dim(C const & A) const
{
	int dim = static_cast<int>( std::floor( std::sqrt( double(A.size()) ) + 0.5 ) );
	assert( dim*dim == int(A.size()) );
	return dim;
}

} /* namespace Algorithms */
} /* namespace elephon */
