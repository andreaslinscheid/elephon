/*	This file AtomDisplacement.cpp is part of elephon.
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

#include "AtomDisplacement.h"
#include <string>
#include <cmath>

namespace elephon
{
namespace LatticeStructure
{

AtomDisplacement::AtomDisplacement()
{

};

AtomDisplacement::AtomDisplacement(
		std::string kind,
		double magnitude,
		std::vector<double> position,
		std::vector<double> direction,
		double gridPrecision,
		bool symmetricDirection)
{
	this->initialize( std::move(kind),magnitude,
			std::move(position),std::move(direction),gridPrecision,symmetricDirection);
};

void AtomDisplacement::initialize(
		std::string kind,
		double magnitude,
		std::vector<double> position,
		std::vector<double> direction,
		double gridPrecision,
		bool symmetricDirection)
{
	kind_ = std::move(kind);
	magnitude_ = magnitude;
	position_=std::move(position);
	for ( auto &xi : position_ )
		xi -= std::floor(xi + 0.5);
	direction_=std::move(direction);
	equivalencePrc_ = gridPrecision;
	treatDirectionSymmetric_ = symmetricDirection;
	this->normalize_direction();
}

void
AtomDisplacement::scale_position(double scaleX, double scaleY, double scaleZ)
{
	assert( position_.size() == 3 );
	position_[0] *= scaleX;
	position_[1] *= scaleY;
	position_[2] *= scaleZ;
}

void AtomDisplacement::transform(Symmetry::SymmetryOperation const& sop)
{
	sop.apply( position_ );
	this->transform_direction(sop);
}

void
AtomDisplacement::transform_direction(Symmetry::SymmetryOperation const& sop)
{
	sop.rotate_cart( direction_ );
	this->normalize_direction(); //symmetry operations in direct coordinates are not necessary unitary
}

double AtomDisplacement::get_prec() const
{
	return equivalencePrc_;
}

std::vector<double> const &
AtomDisplacement::get_position() const
{
	return position_;
}

std::vector<double> const &
AtomDisplacement::get_direction() const
{
	return direction_;
}

std::string AtomDisplacement::get_kind() const
{
	return kind_;
}

double
AtomDisplacement::get_magnitude() const
{
	return magnitude_;
}

bool
AtomDisplacement::is_plus_minus_displ_equivalent() const
{
	return treatDirectionSymmetric_;
}

bool operator< (AtomDisplacement const& d1, AtomDisplacement const& d2)
{
	assert( d1.equivalencePrc_ == d2.equivalencePrc_ );
	assert( d1.treatDirectionSymmetric_ == d2.treatDirectionSymmetric_ );

	//compare names
	int c = d1.kind_.compare( d2.kind_); //zero if equal
	if ( c < 0 )
		return true;
	if ( c > 0 )
		return false;

	auto cmp = [&] (double a, double b ) {
		if ( std::fabs(a-b) < d1.equivalencePrc_ )
			return 0;
		return (a < b ? 1 : -1);
	};

	//Compare position
	for ( int i = 0 ; i < 3; ++i)
		if ( cmp( d1.position_[i],d2.position_[i]) )
			return d1.position_[i] < d2.position_[i];

	//Position is (almost) equal - compare direction
	if ( d1.treatDirectionSymmetric_)
	{
		for ( int i = 0 ; i < 3; ++i)
			if ( cmp( std::abs(d1.direction_[i]),std::abs(d2.direction_[i])) )
				return std::abs(d1.direction_[i]) < std::abs(d2.direction_[i]);
	}
	else
	{
		for ( int i = 0 ; i < 3; ++i)
			if ( cmp( d1.direction_[i],d2.direction_[i]) )
				return d1.direction_[i] < d2.direction_[i];
	}

	//also equal - thus d1<d2 is false
	return false;
}

bool
operator== (AtomDisplacement const& d1, AtomDisplacement const& d2)
{
	return (not (d1 < d2) ) and (not (d2 < d1));
}

void
AtomDisplacement::normalize_direction()
{
	double norm = std::sqrt(direction_[0]*direction_[0] + direction_[1]*direction_[1] + direction_[2]*direction_[2]);
	for ( auto &d : direction_ )
		d /= norm;
}

} /* namespace LatticeStructure */
} /* namespace elephon */
