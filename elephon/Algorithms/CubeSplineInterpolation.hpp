/*	This file CubeSplineInterpolation.hpp is part of elephon.
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
 *  Created on: Jan 7, 2018
 *      Author: A. Linscheid
 */

#include "Algorithms/CubeSplineInterpolation.h"
#include "Auxillary/AlignedVector.h"
#include "Algorithms/LinearAlgebraInterface.h"
#include "Auxillary/CubicPolynomial.h"
#include <cassert>
#include <algorithm>

namespace elephon
{
namespace Algorithms
{

template<typename TMesh,
		typename TData>
CubeSplineInterpolation<TMesh, TData>::CubeSplineInterpolation()
{
	this->clear();
}

template<typename TMesh,
		typename TData>
void
CubeSplineInterpolation<TMesh, TData>::clear()
{
	polys_.clear();
	rangeOfDefinitionBegin = 0.0;
	rangeOfDefinitionEnd = 0.0;
}


template<typename TMesh,
		typename TData>
TData
CubeSplineInterpolation<TMesh, TData>::operator() (
		TMesh x) const
{
	int poly_index = this->find_polynomial_index_in_range(x);
	assert(poly_index < polys_.size());
	TData result;
	polys_[poly_index].evaluate(x, result);
	return result;
}

template<typename TMesh,
		typename TData>
int
CubeSplineInterpolation<TMesh, TData>::find_polynomial_index_in_range(
		TMesh x ) const
{
	assert(polys_.size()>0);
	if ( polys_[0].get_range_min() > x )
		return 0;

	if ( polys_.rbegin()->get_range_max() <= x )
		return static_cast<int>(polys_.size())-1;

	auto it = std::lower_bound(polys_.begin(), polys_.end(), x,
	        [](Auxillary::CubicPolynomial<TMesh,TData> const & lhs, TMesh rhs) -> bool
			  { return lhs.get_range_max() < rhs; });
	// now it points to the first polynom p for which not p.get_range_max() < x

	return std::distance(polys_.begin(), it);
}

template<typename TMesh,
		typename TData>
TData
CubeSplineInterpolation<TMesh, TData>::prime(
		TMesh x) const
{
	int poly_index = this->find_polynomial_index_in_range(x);
	assert(poly_index < polys_.size());
	return polys_[poly_index].evaluate_derivative(x);
}

template<typename TMesh,
		typename TData>
TMesh
CubeSplineInterpolation<TMesh, TData>::min_range() const
{
	return this->rangeOfDefinitionBegin;
}

template<typename TMesh,
		typename TData>
TMesh
CubeSplineInterpolation<TMesh, TData>::max_range() const
{
	return this->rangeOfDefinitionEnd;
}

template<typename TMesh,
		typename TData>
template<class RandomAccessIteratorData, class RandomAccessIteratorMesh>
void
CubeSplineInterpolation<TMesh, TData>::initialize(
		RandomAccessIteratorData dataBegin,
		RandomAccessIteratorData dataEnd,
		RandomAccessIteratorMesh meshBegin,
		RandomAccessIteratorMesh meshEnd,
		TMesh intervalBegin,
		TMesh intervalEnd,
		std::shared_ptr<Auxillary::alignedvector::aligned_vector<TMesh>> splineMatrixPtr)
{
	this->clear();
	assert(rangeOfDefinitionBegin<=rangeOfDefinitionEnd);
	rangeOfDefinitionBegin = intervalBegin;
	rangeOfDefinitionEnd = intervalEnd;

	using std::pow;
	int splineMatrixDim = std::distance(meshBegin, meshEnd);
	assert( splineMatrixDim > 1);
	assert( splineMatrixDim == std::distance(dataBegin, dataEnd));

	// Create the inverse of the spline matrix or take it from input
	Auxillary::alignedvector::aligned_vector<TMesh> inverseSplineMatrix;

	if ( splineMatrixPtr )
	{
		if ( not splineMatrixPtr->empty() )
		{
			assert(std::pow(splineMatrixDim,2) == splineMatrixPtr->size());
			inverseSplineMatrix = std::move(*splineMatrixPtr.get());
		}
		this->compute_inverse_spline_matrix(meshBegin, meshEnd, inverseSplineMatrix);
	}
	else
	{
		this->compute_inverse_spline_matrix(meshBegin, meshEnd, inverseSplineMatrix);
	}

	//create the spline vector
	auto mesh_0 = *(meshBegin + 0);
	auto mesh_1 = *(meshBegin + 1);
	auto mesh_N_m_1 = *(meshBegin + splineMatrixDim - 1);
	auto mesh_N_m_2 = *(meshBegin + splineMatrixDim - 2);
	auto data_0 = *(dataBegin + 0);
	auto data_1 = *(dataBegin + 1);
	auto data_N_m_1 = *(dataBegin + splineMatrixDim - 1);
	auto data_N_m_2 = *(dataBegin + splineMatrixDim - 2);
	Auxillary::alignedvector::aligned_vector<TData> splineVector(splineMatrixDim,0.0);
	splineVector[0] = TData(3.0)*(data_1-data_0)*TData(1.0/pow(mesh_1 - mesh_0 ,2));
	for (int i=1;i<splineMatrixDim-1; ++i)
	{
		auto data_i = *(dataBegin + i);
		auto data_i_m_1 = *(dataBegin + i - 1);
		auto data_i_p_1 = *(dataBegin + i + 1);
		auto mesh_i = *(meshBegin + i);
		auto mesh_i_m_1 = *(meshBegin + i - 1);
		auto mesh_i_p_1 = *(meshBegin + i + 1);
		splineVector[i] = TData(3.0)*((data_i-data_i_m_1)*TData(1.0/pow(mesh_i-mesh_i_m_1,2))
				+(data_i_p_1-data_i)*TData(1.0/pow(mesh_i_p_1-mesh_i,2)));
	}
	splineVector[splineMatrixDim-1] = TData(3.0)*(data_N_m_1-data_N_m_2)
		*TData(1.0/pow(mesh_N_m_1-mesh_N_m_2,2));

	std::vector<TData> vectorOfDerivatives(splineMatrixDim,0.0);
	for (int i=0;i<splineMatrixDim;i++){
		for (int j=0;j<splineMatrixDim;j++)
			vectorOfDerivatives[i] += splineVector[j]*TData(inverseSplineMatrix[i*splineMatrixDim + j]);
	}

	polys_.reserve(splineMatrixDim - 1);
	for (int ipol=0;ipol< splineMatrixDim - 1 ;ipol++){

		// construct a polynomial and insert into the container
		auto data_i = *(dataBegin + ipol);
		auto data_i_p_1 = *(dataBegin + ipol + 1);
		auto mesh_i = *(meshBegin + ipol);
		auto mesh_i_p_1 = *(meshBegin + ipol + 1);
		Auxillary::CubicPolynomial<TMesh,TData> p(
				mesh_i,mesh_i_p_1,data_i,data_i_p_1,
				vectorOfDerivatives[ipol],vectorOfDerivatives[ipol+1]);
		polys_.push_back(std::move(p));
	}

	// if the pointer is set, we place the inverse spline matrix behind it for external further use
	if ( splineMatrixPtr )
		splineMatrixPtr = std::make_shared<Auxillary::alignedvector::aligned_vector<TMesh>>(
				std::move(inverseSplineMatrix));

}

template<typename TMesh,
		typename TData>
template<class VTMesh, class VTData>
void
CubeSplineInterpolation<TMesh, TData>::initialize(
		VTMesh const & mesh,
		VTData const & data)
{
	auto minMaxRange = std::minmax_element(mesh.begin(), mesh.end());
	this->initialize(data.begin(), data.end(),
			mesh.begin(), mesh.end(),
			*minMaxRange.first, *minMaxRange.second);
}

template<typename TMesh,
		typename TData>
template<class RandomAccessIteratorMesh>
void
CubeSplineInterpolation<TMesh, TData>::compute_inverse_spline_matrix(
		RandomAccessIteratorMesh meshBegin,
		RandomAccessIteratorMesh meshEnd,
		Auxillary::alignedvector::aligned_vector<TMesh> & inverseSplineMatrix ) const
{
	int splineMatrixDim = std::distance(meshBegin, meshEnd);
	assert( splineMatrixDim > 1);
	auto mesh_0 = *(meshBegin + 0);
	auto mesh_1 = *(meshBegin + 1);
	auto mesh_N_m_1 = *(meshBegin + splineMatrixDim - 1);
	auto mesh_N_m_2 = *(meshBegin + splineMatrixDim - 2);

	Auxillary::alignedvector::aligned_vector<TMesh>
		splineMatrix(splineMatrixDim*splineMatrixDim,0.0);
	for (int i=1; i<splineMatrixDim-1; ++i)
	{
		auto mesh_i = *(meshBegin + i);
		auto mesh_i_m_1 = *(meshBegin + i - 1);
		auto mesh_i_p_1 = *(meshBegin + i + 1);
		//element i,i-1
		splineMatrix[i*splineMatrixDim+i-1] = 1.0/(mesh_i-mesh_i_m_1);
		//element i,i
		splineMatrix[i*splineMatrixDim+i] = 2.0*(
				1.0/(mesh_i-mesh_i_m_1) + 1.0/(mesh_i_p_1-mesh_i) );
		//element i,+i
		splineMatrix[i*splineMatrixDim+i+1] = 1.0/(mesh_i_p_1-mesh_i);
	}
	// element 0,0
	splineMatrix[0] = 2.0/(mesh_1 - mesh_0 );
	// element 0,1
	splineMatrix[1] = 1.0/(mesh_1 - mesh_0 );
	// element N,N-1
	splineMatrix[(splineMatrixDim-1)*splineMatrixDim+splineMatrixDim-2]
				  = 1.0/(mesh_N_m_1-mesh_N_m_2);
	// element N,N
	splineMatrix[(splineMatrixDim-1)*splineMatrixDim+splineMatrixDim-1]
				  = 2.0/(mesh_N_m_1-mesh_N_m_2);

	// invert the matrix
	Algorithms::LinearAlgebraInterface linalg;
	linalg.inverse(std::move(splineMatrix), inverseSplineMatrix);
}

} /* namespace Algorithms */
} /* namespace elephon */
