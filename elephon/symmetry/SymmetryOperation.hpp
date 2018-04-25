/*	This file SymmetryOperation.hpp is part of elephon.
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
 *  Created on: Jan 16, 2018
 *      Author: A. Linscheid
 */

#include "symmetry/SymmetryOperation.h"
#include "Auxillary/memory_layout_functions.hpp"
#include "Algorithms/LinearAlgebraInterface.h"
#include "Algorithms/helperfunctions.hpp"
#include <cassert>

namespace elephon
{
namespace symmetry
{

template<class VT>
void
SymmetryOperation::rotate_matrix_cart(VT & m) const
{
	typedef typename VT::value_type T;
	assert(m.size()==9);
	auto buffer = m;
	std::fill(buffer.begin(), buffer.end(), T(0));
	for ( int i = 0; i < 3; ++i)
		for ( int j = 0; j < 3; ++j)
			for ( int k = 0; k < 3; ++k)
				for ( int l = 0; l < 3; ++l)
					buffer[i*3+j] += ptgCart[i*3+k]*m[k*3+l]*ptgCart[j*3+l];
}

inline double
SymmetryOperation::get_carth_rot_matrix(int i, int j) const
{
	assert((i>=0)&&(i<3));
	assert((j>=0)&&(j<3));
	return ptgCart[i*3+j];
}

inline double
SymmetryOperation::get_carth_frac_trans(int i) const
{
	assert((i>=0)&&(i<3));
	return fracTransCart[i];
}

inline int
SymmetryOperation::get_lat_rot_matrix (int i, int j) const
{
	assert((i>=0)&&(i<3));
	assert((j>=0)&&(j<3));
	return ptgroup[i*3+j];
}

inline int
SymmetryOperation::get_lat_frac_trans(int i) const
{
	assert((i>=0)&&(i<3));
	return fracTrans[i];
}

template<class VT>
void
SymmetryOperation::rotate_cart( VT & v ) const
{
	assert(v.size()==3);
	this->rotate_cart(v.begin(), v.end());
};

template<typename iterator>
void
SymmetryOperation::rotate_cart( iterator begin, iterator end ) const
{
	assert(std::distance(begin,end)%3 == 0);
	const int nEle = std::distance(begin,end)/3;
	typedef typename std::iterator_traits<iterator>::value_type T;

	for (int iEle = 0 ; iEle < nEle; ++iEle)
	{
		T b[] = {*begin, *(begin+1), *(begin+2)};
		for ( int i = 0; i < 3; ++i)
		{
			*begin = ptgCart[i*3+0]*b[0]+ptgCart[i*3+1]*b[1]+ptgCart[i*3+2]*b[2];
			++begin;
		}
	}
};

template<class Grid>
void
SymmetryOperation::compute_real_space_map(Grid const & grid,
		std::vector<int> & inverseSymOpMap) const
{
	// compute the fractional translation in units of grid points
	std::vector<int> fracGridInv = {-Algorithms::helperfunctions::nint(fracTrans[0]*grid.get_grid_dim()[0]),
									-Algorithms::helperfunctions::nint(fracTrans[1]*grid.get_grid_dim()[1]),
									-Algorithms::helperfunctions::nint(fracTrans[2]*grid.get_grid_dim()[2]) };

	// this provides the lambda that serves as the location mapping
	// NOTE: The following transforms the point in the argument. Thus
	//		it must be applied as the inverse, a.k.a the transposed matrix
	std::vector<int> xyz(3), xyzBuff(3);
	auto indexMapLmbda = [&] (int cnsqGridIndex){
		grid.get_reducible_to_xyz(cnsqGridIndex, xyzBuff);

		// apply symmetry operation on the lattice grid
		std::copy(fracGridInv.begin(), fracGridInv.end(), xyz.begin());
		for (int i = 0 ; i < 3 ; ++i)
			for (int j = 0 ; j < 3 ; ++j)
				xyz[i] += this->get_lat_rot_matrix(i,j)*xyzBuff[j];

		int index = grid.get_xyz_to_reducible_periodic(xyz);
		return index;
	};
	inverseSymOpMap.resize(grid.get_num_points());
	for (int i = 0 ; i < inverseSymOpMap.size(); ++i)
	{
		inverseSymOpMap[i] = indexMapLmbda(i);
	}
}

template<class VT, class Grid>
void
SymmetryOperation::transform_scalar_field_regular_grid(
		Grid const & grid,
		VT & functionData) const
{
	assert(grid.get_num_points() == functionData.size());
	std::vector<int> indexMapRealSpace;
	this->compute_real_space_map(grid, indexMapRealSpace);
	Algorithms::helperfunctions::transform_scalar_field_cart(
			functionData.data(),
			functionData.data()+functionData.size(),
			indexMapRealSpace);
}

template<class VT, class Grid>
void
SymmetryOperation::transform_vector_field_regular_grid(
		Grid const & grid,
		VT & functionData,
		bool transposeFunctionData) const
{
	assert(grid.get_num_points()*3 == functionData.size());
	std::vector<int> indexMapRealSpace;
	this->compute_real_space_map(grid, indexMapRealSpace);
	std::vector<double> ptgCartVec(&ptgCart[0], &ptgCart[0]+sizeof(ptgCart)/sizeof(double));
	if ( not transposeFunctionData )
	{
		Algorithms::helperfunctions::transform_vector_field_cart(
				functionData.data(),
				functionData.data()+functionData.size(),
				ptgCartVec,
				indexMapRealSpace);
	}
	else
	{
		for (int i = 0 ; i < 3 ; ++i)
			Algorithms::helperfunctions::transform_scalar_field_cart(
					functionData.data() + (functionData.size()/3)*i,
					functionData.data() + (functionData.size()/3)*(i+1),
					indexMapRealSpace);

		// now apply the component rotation by block
		Algorithms::LinearAlgebraInterface linalg;
		const int blockSize = (functionData.size()/3);
		auto buffer = functionData;
//		linalg.matrix_matrix_prod(ptgCartVec,buffer,functionData,3,blockSize);
		for (int i = 0 ; i < 3 ; ++i)
		{
			for (int id = 0 ; id < blockSize ; ++id)
				*(functionData.data() + blockSize*i + id) = ptgCart[i*3+0]*(*(buffer.data() + id))
															+ ptgCart[i*3+1]*(*(buffer.data() + blockSize+id))
															+ ptgCart[i*3+2]*(*(buffer.data() + blockSize*2+id));
		}
	}
};

template<class VT>
void
SymmetryOperation::rotate_radial_scalar_data(
		int lMax,
		int numDataPerLM,
		VT & dataToBeTransformed) const
{
	this->rotate_radial_scalar_data(lMax, numDataPerLM, dataToBeTransformed.begin(), dataToBeTransformed.end());
}

template<class interator>
void
SymmetryOperation::rotate_radial_scalar_data(
		int lMax,
		int numDataPerLM,
		interator dataToBeTransformedBegin,
		interator dataToBeTransformedEnd) const
{
	Algorithms::LinearAlgebraInterface linalg;
	assert(radSym_ptr_->size() >= lMax);
	assert(std::distance(dataToBeTransformedBegin, dataToBeTransformedEnd) == std::pow(lMax+1,2)*numDataPerLM);

	typedef typename std::iterator_traits<interator>::value_type T;
	Auxillary::alignedvector::aligned_vector<T> buffer(numDataPerLM*(2*lMax+1));
	for (int l = 0 ; l <= lMax ; ++l)
	{
		int nl = 2*l+1;
		auto const & wignerD = (*radSym_ptr_)[l];
		assert(wignerD.view_as_matrix().size() == nl*nl);
		auto wignerD_ptr = wignerD.view_as_matrix().data();
		auto expansionData_ptr = &( *dataToBeTransformedBegin ) +
				Auxillary::memlayout::angular_momentum_layout(l,-l)*numDataPerLM;
		linalg.call_gemm('t', 'n',	// See declaration comments; the Wigner matrix is transposed because we transform the coefficients, while the
									//		operator acts on the functions.
				nl, numDataPerLM, nl, // W is a (2*l+1) x (2*l+1) matrix, the expansion data we interpret as a (2*l+1) x numDataPerLM matrix.
				decltype(*wignerD_ptr)(1.0), wignerD_ptr, nl,
				expansionData_ptr, numDataPerLM,
				decltype(*wignerD_ptr)(0.0), buffer.data(), numDataPerLM );

		std::copy(buffer.begin(), buffer.begin()+numDataPerLM*nl, expansionData_ptr);
	}
}

template<class interator>
void
SymmetryOperation::rotate_radial_vector_data(
		int lMax,
		int numDataPerLM,
		interator dataToBeTransformedBegin,
		interator dataToBeTransformedEnd) const
{
	const int blockSize = std::pow(lMax+1,2)*numDataPerLM;
	assert(std::distance(dataToBeTransformedBegin, dataToBeTransformedEnd) == 3*blockSize);

	// First we transform the data in each channel, then we apply the symmetry operation to the channels
	this->rotate_radial_scalar_data(lMax,numDataPerLM, dataToBeTransformedBegin, dataToBeTransformedBegin + blockSize);
	this->rotate_radial_scalar_data(lMax,numDataPerLM, dataToBeTransformedBegin + blockSize, dataToBeTransformedBegin + 2*blockSize);
	this->rotate_radial_scalar_data(lMax,numDataPerLM, dataToBeTransformedBegin + 2*blockSize, dataToBeTransformedEnd);

	typedef typename std::iterator_traits<interator>::value_type T;
	Auxillary::alignedvector::aligned_vector<T> buffer(dataToBeTransformedBegin, dataToBeTransformedEnd);
	for (int i = 0 ; i < 3 ; ++i)
	{
		for (int id = 0 ; id < blockSize ; ++id)
			*(dataToBeTransformedBegin + blockSize*i + id) = *(buffer.data() + id)*ptgCart[i*3+0]
														+ *(buffer.data() + blockSize+id)*ptgCart[i*3+1]
														+ *(buffer.data() + blockSize*2+id)*ptgCart[i*3+2];
	}
}

} /* namespace symmetry */
} /* namespace elephon */
