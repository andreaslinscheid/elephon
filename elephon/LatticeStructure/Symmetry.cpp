/*	This file Symmetry.cpp is part of elephon.
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
 *  Created on: May 15, 2017
 *      Author: A. Linscheid
 */

#include "Symmetry.h"
#include <assert.h>
#include <set>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <map>
#include <algorithm>

namespace elephon
{
namespace LatticeStructure
{

Symmetry::Symmetry()
{
	//This sets the identity symmetry group
	symmetries_ = {	1 , 0 , 0,
					0 , 1 , 0,
					0 , 0 , 1};
	fractTrans_ = std::vector<double>({0.0, 0.0, 0.0});
	lattice_.initialize( std::vector<double>({1,0,0,0,1,0,0,0,1}) );
	this->initialize(1e-6,symmetries_,fractTrans_,lattice_,false);
}

void
Symmetry::initialize(
		double symmPrec,
		std::vector<int> symmetries,
		std::vector<double> fractionalTranslations,
		LatticeStructure::LatticeModule lattice,
		bool hasTimeReversal)
{
	assert( symmetries.size() % 9 == 0 );
	assert( fractionalTranslations.size() % 3 == 0 );
	assert( (fractionalTranslations.size() / 3) == (symmetries.size() / 9) );
	symmetries_ = std::move( symmetries );
	numSymmetries_ = symmetries_.size()/9;
	fractTrans_ = std::move(fractionalTranslations);
	symmPrec_ = symmPrec;
	lattice_ = std::move(lattice);
	hasTimeReversal_ = hasTimeReversal;

	//This transforms S' = B S A into cartesian space
	// as S = A S' B
	auto A = lattice_.get_latticeMatrix();
	auto B = lattice_.get_reciprocal_latticeMatrix();

	auto matMult = [] (double const * A, double const * B, double * result){
		std::fill(result,result+9,0.0);
		for ( int i = 0 ; i < 3 ; ++i)
			for ( int j = 0 ; j < 3 ; ++j)
				for ( int k = 0 ; k < 3 ; ++k)
					result[i*3+j] += A[i*3+k]*B[k*3+j];
	};

	auto transform = [&A,&B] (int const * symOp, double * cartSymOp ) {
		std::fill(cartSymOp,cartSymOp+9, 0.0 );
		for ( int i = 0 ; i < 3 ; ++i)
			for ( int j = 0 ; j < 3 ; ++j)
				for ( int k = 0 ; k < 3 ; ++k)
					for ( int l = 0 ; l < 3 ; ++l)
						cartSymOp[i*3+j] += A[i*3+k]*symOp[k*3+l]*B[l*3+j];
	};

	for (int isym = 0 ; isym < numSymmetries_; ++isym)
		if ( std::vector<int>( &symmetries_[9*isym], &symmetries_[9*isym]+9)
				== std::vector<int>({1,0,0,0,1,0,0,0,1}) )
		{
			idIndex_ = isym;
			break;
		}

	hasInversion_ = false;
	for (int isym = 0 ; isym < numSymmetries_; ++isym)
		if ( std::vector<int>( &symmetries_[9*isym], &symmetries_[9*isym]+9)
				== std::vector<int>({-1,0,0, 0,-1,0, 0,0,-1}) )
		{
			hasInversion_ = true;
			break;
		}

	//In case we don't have inversion but time reversal the number of symmetry operations in
	//reciprocal space is different.
	if ( (not hasInversion_) and hasTimeReversal_ )
		numRotations_ = 2*numSymmetries_;
	else
		numRotations_ = numSymmetries_;

	//For fast lookups, we define an ordering according to the function below.
	//It takes two symmetry operations 'a' and 'b' and decides which one is smaller by an element wise comparison
	//of the nearest integer values. We add an additional element to the front, before
	//the elements of the rotation that flags time reversal symmetry
	auto cmpSyms = [&] ( std::vector<double> const & a, std::vector<double> const & b )
	{
		assert( (a.size()==10) && (b.size()==10) );
		for ( int i = 0 ; i < 10; ++i)
			if ( std::abs(a[i]-b[i]) > symmPrec_ )
				return a[i]<b[i];
		return false;
	};
	std::map< std::vector<double>, int, decltype(cmpSyms)> symmetryCartesianSet(cmpSyms);
	std::vector<double> s_cart(10);
	std::vector<double> inversion = {-1,0,0, 0,-1,0, 0,0,-1};
	symmetriesCartesian_.resize(numRotations_*9);
	for (int isym = 0 ; isym < numSymmetries_; ++isym)
	{
		s_cart[0] = -1.0; // < 0 flags normal symmetry operation
		transform(&symmetries_[9*isym],&symmetriesCartesian_[9*isym]);
		std::copy(&symmetriesCartesian_[9*isym],&symmetriesCartesian_[9*isym]+9,&s_cart[1]);
		symmetryCartesianSet.insert( std::make_pair(s_cart, isym ) );

		//In case of time reversal symmetry but no inversion we add the generated additional
		//rotation into the set. They get a flag which orders them to end of the map.
		if ( (not hasInversion_) and hasTimeReversal_ )
		{
			s_cart[0] = 1.0; // > 0 flags time reversal enabled symmetry operation
			transform(&symmetries_[9*isym],&s_cart[1]);
			matMult(&s_cart[1],inversion.data(),&symmetriesCartesian_[9*(isym+numSymmetries_)]);
			std::copy(&symmetriesCartesian_[9*(isym+numSymmetries_)],
					&symmetriesCartesian_[9*(isym+numSymmetries_)]+9,&s_cart[1]);
			symmetryCartesianSet.insert( std::make_pair(s_cart, isym ) );
		}

	}
	assert( int(symmetryCartesianSet.size()) == numRotations_ );

	//ToDo This does not take into account fractional translations.
	//		This must be handled with more care because a supercell build of
	//		a unit cell with fractional translations can have identity operations
	//		with fractional translations.
	inverseMap_.resize(numRotations_);
	multiplicationTable_.resize( numRotations_*numRotations_);
	for (auto sym : symmetryCartesianSet)
	{
		//find the group multiplication table
		std::vector<double> prod(10);
		int inverseIndex = -1;
		for (auto secondSym : symmetryCartesianSet)
		{
			matMult(&sym.first[1],&secondSym.first[1],&prod[1]);
			prod[0] = sym.first[0] == secondSym.first[0] ? -1.0 : 1.0; //Do we have an even number of inversions in this product?
			auto ret = symmetryCartesianSet.find( prod );
			if ( ret == symmetryCartesianSet.end() )
				throw std::runtime_error("Symmetries do not form a group!");
			multiplicationTable_[sym.second*numRotations_+secondSym.second] = ret->second;
			if ( ret->second == this->get_identity_index() )
				inverseIndex = secondSym.second;
		}

		if ( inverseIndex < 0 )
			throw std::runtime_error("Found no inverse matrix");
		inverseMap_[sym.second] = inverseIndex;
	}

	//Append the symmetry operations and fractional translations
	if ( (not hasInversion_) and hasTimeReversal_ )
	{
		symmetries_.reserve(numRotations_*9);
		fractTrans_.reserve(numRotations_*3);
		for (int isym = 0; isym < numSymmetries_; ++isym)
		{
			// S' = -I * S, so we multiply each element by -1
			for (int i = 0; i < 9; ++i)
				symmetries_.push_back(-symmetries_[isym*9+i]);
			for (int i = 0; i < 3; ++i)
				fractTrans_.push_back(fractTrans_[isym*3+i]);
		}
	}

	isSymmorphic_ = std::vector<bool>( numRotations_, true );
	for ( int irot = 0 ; irot < numRotations_ ; ++irot  )
		isSymmorphic_[irot] = (std::abs(fractTrans_[irot*3]) < symmPrec_ )
							and (std::abs(fractTrans_[irot*3+1]) < symmPrec_ )
							and (std::abs(fractTrans_[irot*3+2]) < symmPrec_ ) ;

	fractTransCartesian_.resize( numRotations_*3 );
	for (int isym = 0; isym < numRotations_; ++isym)
		for (int i = 0; i < 3; ++i)
			fractTransCartesian_[isym*3+i]
					= A[i*3+0]*fractTrans_[isym*3+0] + A[i*3+1]*fractTrans_[isym*3+1] + A[i*3+2]*fractTrans_[isym*3+2];

	//TODO We should make sure that the lattice and the symmetries are compatible
}

void
Symmetry::set_reciprocal_space_sym()
{
	if ( isReciprocalSpace_ )
		return;
	fractTransStore_ = fractTrans_;
	std::fill(fractTrans_.begin(),fractTrans_.end(),0.0);
	std::fill(fractTransCartesian_.begin(),fractTransCartesian_.end(),0.0);
	isReciprocalSpace_ = true;
}

void
Symmetry::apply(int isym, std::vector<double> & field, bool latticePeriodic) const
{
	this->apply(isym,field.begin(),field.end(),latticePeriodic);
}

void
Symmetry::apply(int isym, std::vector<double>::iterator fieldBegin,
					std::vector<double>::iterator fieldEnd,
					bool latticePeriodic) const
{
	assert( std::distance(fieldBegin,fieldEnd)%3 == 0 );
	assert( isym < numRotations_ );

	std::vector<double> buff(3);

	int numComponents = std::distance(fieldBegin,fieldEnd)/3;
	for ( int ic = 0; ic < numComponents; ++ic)
	{
		std::copy(fieldBegin+ic*3,fieldBegin+(ic+1)*3,std::begin(buff));
		if ( not isReciprocalSpace_ )
			for ( int xi = 0; xi < 3; ++xi)
			{
				*(fieldBegin+ic*3+xi) = symmetries_[(isym*3+xi)*3+0]*buff[0]
								+symmetries_[(isym*3+xi)*3+1]*buff[1]
								+symmetries_[(isym*3+xi)*3+2]*buff[2]
								+fractTrans_[isym*3+xi];
				if ( latticePeriodic )
					*(fieldBegin+ic*3+xi) -= std::floor( *(fieldBegin+ic*3+xi)+0.5);
			}
		else
			for ( int xi = 0; xi < 3; ++xi)
			{
				*(fieldBegin+ic*3+xi) = symmetries_[(isym*3+0)*3+xi]*buff[0]
								+symmetries_[(isym*3+1)*3+xi]*buff[1]
								+symmetries_[(isym*3+2)*3+xi]*buff[2];
				if ( latticePeriodic )
					*(fieldBegin+ic*3+xi) -= std::floor( *(fieldBegin+ic*3+xi)+0.5);
			}
	}
}

void
Symmetry::apply_cartesian(int isym, std::vector<double>::iterator fieldCartBegin,
					std::vector<double>::iterator fieldCartEnd) const
{
	assert( std::distance(fieldCartBegin,fieldCartEnd)%3 == 0 );
	assert( isym < numRotations_ );

	if ( isReciprocalSpace_ )
		//in cartesian space, symmetries are unitary
		isym = this->get_index_inverse(isym);

	std::vector<double> buff(3);

	int numComponents = std::distance(fieldCartBegin,fieldCartEnd)/3;
	for ( int ic = 0; ic < numComponents; ++ic)
	{
		std::copy(fieldCartBegin+ic*3,fieldCartBegin+(ic+1)*3,std::begin(buff));
		for ( int xi = 0; xi < 3; ++xi)
		{
			*(fieldCartBegin+ic*3+xi) = symmetriesCartesian_[(isym*3+xi)*3+0]*buff[0]
							+symmetriesCartesian_[(isym*3+xi)*3+1]*buff[1]
							+symmetriesCartesian_[(isym*3+xi)*3+2]*buff[2]
							+fractTransCartesian_[isym*3+xi];
		}
	}
}

void
Symmetry::rotate_cartesian(int isym, std::vector<double>::iterator fieldCartBegin,
					std::vector<double>::iterator fieldCartEnd) const
{
	assert( std::distance(fieldCartBegin,fieldCartEnd)%3 == 0 );
	assert( isym < numRotations_ );

	if ( isReciprocalSpace_ )
		//in cartesian space, symmetries are unitary
		isym = this->get_index_inverse(isym);

	std::vector<double> buff(3);

	int numComponents = std::distance(fieldCartBegin,fieldCartEnd)/3;
	for ( int ic = 0; ic < numComponents; ++ic)
	{
		std::copy(fieldCartBegin+ic*3,fieldCartBegin+(ic+1)*3,std::begin(buff));
		for ( int xi = 0; xi < 3; ++xi)
		{
			*(fieldCartBegin+ic*3+xi) = symmetriesCartesian_[(isym*3+xi)*3+0]*buff[0]
							+symmetriesCartesian_[(isym*3+xi)*3+1]*buff[1]
							+symmetriesCartesian_[(isym*3+xi)*3+2]*buff[2];
		}
	}
}

void
Symmetry::rotate_matrix_cartesian(int isym,
		std::vector<double>::iterator matrixFieldCartBegin,
		std::vector<double>::iterator matrixFieldCartEnd) const
{
	int elem = std::distance(matrixFieldCartBegin,matrixFieldCartEnd);
	assert( elem%9 == 0 );
	assert( isym < numRotations_ );

	if ( isReciprocalSpace_ )
		//in cartesian space, symmetries are unitary
		isym = this->get_index_inverse(isym);

	std::vector<double> b(9);

	int numComponents = elem/9;
	for ( int ic = 0; ic < numComponents; ++ic)
	{
		std::copy(matrixFieldCartBegin+ic*9,matrixFieldCartBegin+(ic+1)*9,std::begin(b));
		std::fill(matrixFieldCartBegin+ic*9,matrixFieldCartBegin+(ic+1)*9, 0.0 );
		auto G = &symmetriesCartesian_[isym*9];
		for ( int i = 0; i < 3; ++i)
			for ( int j = 0; j < 3; ++j)
				for ( int k = 0; k < 3; ++k)
					for ( int l = 0; l < 3; ++l)
						*(matrixFieldCartBegin+ic*9+i*3+l) += G[i*3+j]*b[j*3+k]*G[l*3+k];
	}
}

double
Symmetry::get_symmetry_prec() const
{
	return symmPrec_;
}

int
Symmetry::get_num_symmetries() const
{
	return isReciprocalSpace_ ? numRotations_ : numSymmetries_ ;
}

void
Symmetry::small_group(std::vector<double> const& point)
{
	//apply all symmetry operations and discard those who map point somewhere else
	std::vector<int> dropIndices;
 	for (int isym = 0 ; isym < numRotations_; ++isym)
	{
		auto rotcpy = point;
		this->apply(isym, rotcpy.begin(), rotcpy.end(), false );
		for ( int i = 0; i < int(rotcpy.size()); ++i)
			if ( std::abs(rotcpy[i]-point[i]) > symmPrec_ )
			{
				dropIndices.push_back(isym);
				break;
			}
	}
 	this->symmetry_reduction(dropIndices);
}

void
Symmetry::small_group_cart(std::vector<double> const& pointCartCoords)
{
	//apply all symmetry operations and discard those who map point somewhere else
	std::vector<int> dropIndices;
 	for (int isym = 0 ; isym < numRotations_; ++isym)
	{
		auto rotcpy = pointCartCoords;
		this->apply_cartesian(isym, rotcpy.begin(), rotcpy.end() );
		for ( int i = 0; i < int(rotcpy.size()); ++i)
			if ( std::abs(rotcpy[i]-pointCartCoords[i]) > symmPrec_ )
			{
				dropIndices.push_back(isym);
				break;
			}
	}
 	this->symmetry_reduction(dropIndices);
}

bool
Symmetry::is_symmorphic(int isym) const
{
	assert(isym < numRotations_);
	return isSymmorphic_[isym];
}

int
Symmetry::get_index_inverse(int isym) const
{
	return inverseMap_[isym];
}

int
Symmetry::get_group_product(int isym1, int isym2) const
{
	assert( (isym1 < numRotations_) && (isym2 < numRotations_) );
	return multiplicationTable_[isym1*numRotations_+isym2];
}

bool
Symmetry::has_inversion() const
{
	return hasInversion_;
}

bool
Symmetry::is_inversion(int isym) const
{
	assert( isym < numRotations_ );
	return (isym >= numSymmetries_) and ( isym < numRotations_);
}

void
Symmetry::symmetry_reduction( std::vector<int> indicesDropped)
{
	//In case of time reversal expanded symmetry sets, remove those from the list
	//The new initialization will put them back in.
	if ( numRotations_ != numSymmetries_)
	{
		auto it = std::remove_if(indicesDropped.begin(), indicesDropped.end(),
				std::bind2nd(std::greater_equal<int>(), numSymmetries_));
		indicesDropped.erase (it, indicesDropped.end());
	}

	if ( indicesDropped.empty() )
		return;

	//remove possible duplicates and make sure all indices appear
	std::set<int> drop(indicesDropped.begin(),indicesDropped.end());
	assert( ((*drop.rbegin()) < numSymmetries_) && ((*drop.begin()) >= 0) );
	if ( drop.find(idIndex_) !=  drop.end() )
		throw std::logic_error("Requested to remove the identity element from the symmetry ... ");

	int numNewSym = numSymmetries_ - static_cast<int>( drop.size() );
	std::vector<int> newPtGrpSym( numNewSym*9 );
	std::vector<double> newFractSym( numNewSym*3 );
	int nis = 0;
	for ( int isym = 0 ; isym < numSymmetries_; ++isym)
		if ( drop.find(isym) == drop.end() ) // we keep this one
		{
			std::copy(symmetries_.begin()+isym*9,symmetries_.begin()+(isym+1)*9,newPtGrpSym.begin()+nis*9);
			std::copy(fractTrans_.begin()+isym*3,fractTrans_.begin()+(isym+1)*3,newFractSym.begin()+nis*3);
			nis++;
		}
	assert( numNewSym == nis);
	this->initialize( symmPrec_, newPtGrpSym, newFractSym, lattice_ , hasTimeReversal_);
}

Symmetry::SymmetryOperation
Symmetry::get_sym_op( int isym ) const
{
	assert( isym < numRotations_ );
	Symmetry::SymmetryOperation res;
	for ( int i = 0 ; i < 3 ; ++ i)
	{
		res.fracTrans[i] = (isReciprocalSpace_ ? -1.0 : 1.0)*fractTrans_[isym*3+i];
		for ( int j = 0 ; j < 3 ; ++j)
			res.ptgroup[i*3+j] = symmetries_[isym*9 + ((isReciprocalSpace_ ? j*3+i : i*3+j)) ];
	}

	if ( isReciprocalSpace_ )
		//in cartesian space, symmetries are unitary
		isym = this->get_index_inverse(isym);
	std::copy( &symmetriesCartesian_[isym*9],&symmetriesCartesian_[isym*9]+9,res.ptgCart);
	return res;
}


void
Symmetry::Sop::apply( std::vector<double> & v, bool latticePeriodic ) const
{
	assert(v.size()==3);
	auto b = v;
	for ( int i = 0; i < 3; ++i)
	{
		v[i] = ptgroup[i*3+0]*b[0]+ptgroup[i*3+1]*b[1]+ptgroup[i*3+2]*b[2]+fracTrans[i];
		if ( latticePeriodic )
			v[i] -= std::floor(v[i]+0.5);
	}
};

void
Symmetry::Sop::rotate_cart( std::vector<double> & v ) const
{
	assert(v.size()==3);
	auto b = v;
	for ( int i = 0; i < 3; ++i)
		v[i] = ptgCart[i*3+0]*b[0]+ptgCart[i*3+1]*b[1]+ptgCart[i*3+2]*b[2];
};

int
Symmetry::get_identity_index() const
{
	return idIndex_;
}

bool
Symmetry::is_reci() const
{
	return isReciprocalSpace_;
}

std::vector<double>
Symmetry::get_fractional_translation(int isym) const
{
	std::vector<double> result(3);
	if ( isReciprocalSpace_ )
	{
		assert( int(fractTransStore_.size()) == numSymmetries_*3);
		std::copy(fractTransStore_.begin()+isym*3,fractTransStore_.begin()+(isym+1)*3,result.begin());
	}
	else
	{
		assert( int(fractTrans_.size()) == numSymmetries_*3);
		std::copy(fractTrans_.begin()+isym*3,fractTrans_.begin()+(isym+1)*3,result.begin());
	}
	return result;
}

LatticeStructure::LatticeModule const &
Symmetry::get_lattice() const
{
	return lattice_;
}

void
Symmetry::reset_lattice(LatticeStructure::LatticeModule lattice)
{
	//simply recalculate everything ...
	std::vector<int> newPtGrpSym(9*numSymmetries_);
	std::vector<double> newFractSym(3*numSymmetries_);
	for ( int isym = 0 ; isym < numSymmetries_; ++isym)
	{
		std::copy(symmetries_.begin()+isym*9,symmetries_.begin()+(isym+1)*9,newPtGrpSym.begin()+isym*9);
		std::copy(fractTrans_.begin()+isym*3,fractTrans_.begin()+(isym+1)*3,newFractSym.begin()+isym*3);
	}
	this->initialize( symmPrec_, newPtGrpSym, newFractSym, lattice, hasTimeReversal_);
}

} /* namespace LatticeStructure */
} /* namespace elephon */