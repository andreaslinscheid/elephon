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
#include <cmath>
#include <assert.h>
#include <stdexcept>

#include <iostream>

namespace elephon
{
namespace Algorithms
{

template<class C>
void
LinearAlgebraInterface::pseudo_inverse(C A, int m, int n,
		C & pinvA, double cutoff)
{
	typedef typename C::value_type T;
	//This routine needs to buffer a few matrices.
	//In a kinda hacky fashion, we put them in the workspace, too.
	typedef typename detail::ComplexTypeTrait<typename C::value_type>::type realT;

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
LinearAlgebraInterface::null_space(std::vector<T> A, int m, int n,
		int & kerDim, std::vector<T> & nullA, double tol)
{
	std::vector<T> U(m*m);
	std::vector<T> VT(n*n);
	std::vector<typename detail::ComplexTypeTrait<T>::type > sv(std::min(n,m));
	this->svd(std::move(A),m,n,U,VT,sv);

	kerDim = 0;
	for ( int i = 0 ; i < sv.size(); ++i)
		if ( sv[i] < tol)
			kerDim++;

	nullA.resize(kerDim*n);
	for ( int i = n-kerDim ; i < n; ++i)
		for ( int j = 0 ; j < n; ++j)
		{
			nullA[(i-n+kerDim)*n+j] = VT[i*n+j];
		}
}

template<class VT1, class VT2,class VT3>
void
LinearAlgebraInterface::matrix_matrix_prod(VT1 const & A,
		VT2 const & B,
		VT3 & ATimesB, int m, int n) const
{
	typedef typename VT1::value_type T;
	typedef typename VT2::value_type T2;
	typedef typename VT3::value_type T3;
	static_assert(std::is_same<T,T2>::value, "Value types for matrix_matrix_prod of first and second argument must match");
	static_assert(std::is_same<T,T3>::value, "Value types for matrix_matrix_prod of first and third argument must match");
	int k = int(A.size())/m;
	assert( k == int(B.size())/n );
	assert( ATimesB.size() == m * n );
	this->call_gemm('n','n',m,n,k,T(1.0),A.data(),k,B.data(),n,T(0.0),ATimesB.data(),n);
}

template<typename VT>
void
LinearAlgebraInterface::matrix_vector_prod(
		VT const & A,
		VT const & B,
		VT & ATimesB) const
{
	typedef typename VT::value_type T;
	int n = static_cast<int>(B.size());
	assert(n > 0);
	int m = static_cast<int>(A.size())/n;
	assert( static_cast<int>(A.size()) == m * n );
	ATimesB.resize( m );
	this->call_gemv('n',m,n,T(1.0),A.data(),n,B.data(),1,T(0.0),ATimesB.data(),1);
}

template<typename T>
void
LinearAlgebraInterface::svd(std::vector<T> A, int m, int n,
		std::vector<T> & U,
		std::vector<T> & VT,
		std::vector< typename detail::ComplexTypeTrait< std::complex<T> >::type > & sv)
{
	typedef typename detail::ComplexTypeTrait<T>::type realT;
	U.resize(m*m);
	VT.resize(n*n);
	sv.resize(std::min(m,n));

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
	int nElemChar = sizeof(T)*lwork;
	workbuffer_.resize( nElemChar );
	T * work = reinterpret_cast< T * >( &workbuffer_[0] );

	if ( sizeof(realT) != sizeof(T) )
		//We are running the complex version
		rWork_.resize( sizeof(realT)*(std::max(1,std::min(m,n)*std::max(5*std::min(m,n)+7,2*std::max(m,n)+2*std::min(m,n)+1))) );
		//The above formula (mess) is from https://software.intel.com/en-us/mkl-developer-reference-fortran-gesdd

	//The complex version needs a parameter rwork
	realT * rwork =  reinterpret_cast< realT * >( rWork_.data() );

	//Left singular vectors (U) are stored column wise, right singular vectors (vT) are stored row wise
	this->call_gesdd( LAPACK_ROW_MAJOR, 'A', m, n, A.data(), lda, sv.data(), U.data(), ldu, VT.data(),
			ldvt, work, rwork, lwork, IPIV_.data()  );
}

template<typename VT>
void
LinearAlgebraInterface::transpose_square_matrix(VT & A) const
{
	int dim = this->square_matrix_dim(A);
	for (int i = 0 ; i < dim ; ++i)
		for (int j = i+1 ; j < dim ; ++j)
		{
			typename VT::value_type tmp = A[i*dim+j];
			A[i*dim+j] = A[j*dim+i];
			A[j*dim+i] = tmp;
		}
}

template<typename VT>
void
LinearAlgebraInterface::conjugate_square_matrix(VT & A) const
{
	this->transpose_square_matrix(A);
	for (auto & a : A )
		a = std::conj(a);
}

template<class C>
void
LinearAlgebraInterface::inverse(C A, C & invA )
{
	typedef typename std::remove_pointer<decltype(A.data())>::type T;
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
	int dim = this->square_matrix_dim( matrix );
	eigenvectors = matrix;
	eigenvalues.resize(dim);

	char v = ( comEV ?  'V' : 'N' );

	char dataLoc = ( data_upper ? 'U' : 'L' );


	// first attempt to diagonalize with heevr
	// according to Intel, heevr should be the default for diagonalizing hermitian matrices
	int info = 0;
	int liwork = -1;
	int lrwork = -1;
	int lwork = -1;
	std::vector<int> & iwork  = IPIV_; // use same naming convention as in the interface
	int iwork_query;
	T rwork_query;
	std::complex<T> work_query;
	int m;
	std::vector<int> isuppz(2*dim);

	char range = 'A'; // we want all eigenvalues
	T vl = 0, vu = 0;
	int il = 0, iu = 0;
	T abstol = 1e-8; // all eigenproblems we are concerned with in this code are O(1)

	// Query optimal working array size
	info = this->call_heevr( LAPACK_ROW_MAJOR, v, dataLoc, range, dim, matrix.data(), dim, vl,
							 vu, il, iu, abstol, &m, eigenvalues.data(), eigenvectors.data(), dim, isuppz.data(),
							 &work_query, lwork, &rwork_query, lrwork,
							 &iwork_query, liwork );

	int optimal_size = static_cast<int>( std::real(work_query) ) ;
	if ( workbuffer_.size() < sizeof(std::complex<T>)*optimal_size )
		workbuffer_.resize( sizeof(std::complex<T>)*optimal_size );

	liwork = (int)iwork_query;
	if ( iwork.size() < sizeof(T)*liwork )
		iwork.resize(sizeof(T)*liwork);

	lrwork = (int)rwork_query;
	if ( rWork_.size() < sizeof(T)*lrwork )
		rWork_.resize(sizeof(T)*lrwork);

	info = this->call_heevr( LAPACK_ROW_MAJOR, v, dataLoc, range, dim, matrix.data(), dim, vl,
							 vu, il, iu, abstol, &m, eigenvalues.data(), eigenvectors.data(), dim, isuppz.data(),
							 reinterpret_cast< std::complex<T>* >(workbuffer_.data()), optimal_size,
							 reinterpret_cast< T* >(rWork_.data()), lrwork,
							 iwork.data(), liwork );

	// backup if heevr fails
	if ( info != 0 )
	{
		int rworksize = ( 3*dim-2 > 1 ? 3 * dim - 2 : 1 ) ;
		if ( rWork_.size() < sizeof(T)*rworksize )
			rWork_.resize(sizeof(T)*rworksize);

		std::complex<T> work_query;
		this->call_heev(
				LAPACK_ROW_MAJOR,
				v, dataLoc,
				dim, eigenvectors.data(), dim,
				eigenvalues.data(),&work_query,-1,
				reinterpret_cast< T* >(rWork_.data()));
		size_t optimal_size = static_cast<size_t>( std::real(work_query) ) ;
		if ( workbuffer_.size() < sizeof(std::complex<T>)*optimal_size )
			workbuffer_.resize( sizeof(std::complex<T>)*optimal_size );

		int lwork = static_cast<int>(optimal_size);
		info = this->call_heev(
				LAPACK_ROW_MAJOR,
				v, dataLoc,
				dim, eigenvectors.data(), dim,
				eigenvalues.data(),
				reinterpret_cast< std::complex<T>* >(workbuffer_.data()),
				lwork,
				reinterpret_cast< T* >(rWork_.data()));

		this->check_library_info(info, "heev diag");
	}

	//Ensure phase convention on eigenvectors (real part of the last component is non-negative)
	for ( int ib = 0 ; ib < dim ; ++ib )
		if ( eigenvectors[(dim-1)*dim + ib].real() < 0 )
			for ( int ibp = 0 ; ibp < dim ; ++ibp )
				eigenvectors[ibp*dim + ib] *= std::complex<T>(-1.0);
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
