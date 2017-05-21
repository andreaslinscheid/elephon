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

namespace elephon
{
namespace LatticeStructure
{

Symmetry::Symmetry()
{
	//This sets the identity symmetry
	symmetries_ = {	1 , 0 , 0,
					0 , 1 , 0,
					0 , 0 , 1};
	numSymmetries_ = 1;
	fractTrans_ = {0.0, 0.0, 0.0};
	symmPrec_ = 1e-6;
}

void
Symmetry::initialize(
		double symmPrec,
		std::vector<int> symmetries,
		std::vector<double> fractionalTranslations )
{
	assert( symmetries.size() % 9 == 0 );
	assert( fractionalTranslations.size() % 3 == 0 );
	symmetries_ = std::move( symmetries );
	numSymmetries_ = symmetries_.size()/9;
	fractTrans_ = std::move(fractionalTranslations);
	symmPrec_ = symmPrec;

	//TODO in principle here we should check if the input is sensible, e.g. if the symmetries
	//form a group.
	//ToDo make sure the identity is in place index 0
}

void
Symmetry::set_reciprocal_space_sym()
{
	std::fill(fractTrans_.begin(),fractTrans_.end(),0.0);
	isReciprocalSpace_ = true;
}

void
Symmetry::apply(int isym, std::vector<double> & field, bool latticePeriodic) const
{
	assert( field.size()%3 == 0 );
	assert( isym < numSymmetries_ );

	std::vector<double> buff(3);

	int numComponents = static_cast<int>(field.size())/3;
	for ( int ic = 0; ic < numComponents; ++ic)
	{
		std::copy(std::begin(field)+ic*3,std::begin(field)+(ic+1)*3,std::begin(buff));
		for ( int xi = 0; xi < 3; ++xi)
		{
			field[ic*3+xi] = symmetries_[(isym*3+xi)*3+0]*buff[0]
							+symmetries_[(isym*3+xi)*3+1]*buff[1]
							+symmetries_[(isym*3+xi)*3+2]*buff[2]
							+fractTrans_[isym*3+xi];
			if ( latticePeriodic )
				field[ic*3+xi] -= std::floor(field[ic*3+xi]);
		}
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
	return numSymmetries_;
}

void
Symmetry::symmetry_reduction( std::vector<int> const& indicesDropped)
{
	//remove possible duplicates and make sure all indices appear
	std::set<int> drop(indicesDropped.begin(),indicesDropped.end());
	assert( ((*drop.rbegin()) < numSymmetries_) && ((*drop.begin()) >= 0) );

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
	this->initialize( symmPrec_, newPtGrpSym, newFractSym );
}

Symmetry::SymmetryOperation
Symmetry::get_sym_op( int isym ) const
{
	assert( isym < numSymmetries_ );
	Symmetry::SymmetryOperation res;
	std::copy( &fractTrans_[isym*3],&fractTrans_[isym*3]+3,res.fracTrans);
	std::copy( &symmetries_[isym*9],&symmetries_[isym*9]+9,res.ptgroup);
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

int
Symmetry::get_identity_index() const
{
	return 0;
}

bool
Symmetry::is_reci() const
{
	return isReciprocalSpace_;
}

} /* namespace LatticeStructure */
} /* namespace elephon */
