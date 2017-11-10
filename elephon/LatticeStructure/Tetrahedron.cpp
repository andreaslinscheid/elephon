/*	This file Tetrahedron.cpp is part of elephon.
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
 *  Created on: Nov 2, 2017
 *      Author: A. Linscheid
 */

#include "LatticeStructure/Tetrahedron.h"

namespace elephon
{
namespace LatticeStructure
{

Tetrahedron::Tetrahedron(
		std::vector<int> cornerIndicesData,
		std::vector<int> cornerIndicesExtended,
		std::vector<int> cornerIndicesExtendedReducible,
		std::shared_ptr<const ExtendedSymmetricGrid> extendedGrid)
{
	assert(cornerIndicesData.size() == 4);
	assert(cornerIndicesExtended.size() == 4);
	assert(cornerIndicesExtendedReducible.size() == 4);

	std::multimap<int,int> sorter;
	for ( int i = 0 ; i < 4 ; ++i)
		sorter.insert(std::make_pair(cornerIndicesExtended[i],i));
	cornerIndicesExtended_.clear();
	cornerIndicesExtendedReducible_.clear();
	cornerIndicesData_.clear();
	cornerIndicesExtended_.reserve(4);
	cornerIndicesExtendedReducible_.reserve(4);
	cornerIndicesData_.reserve(4);
	for ( auto s : sorter)
	{
		cornerIndicesExtended_.push_back( s.first );
		cornerIndicesExtendedReducible_.push_back(cornerIndicesExtendedReducible[s.second]);
		cornerIndicesData_.push_back(cornerIndicesData[s.second]);
	}

	extendedGrid_ = extendedGrid;
}

int
Tetrahedron::get_multiplicity() const
{
	return multiplicity_;
}

void
Tetrahedron::set_multiplicity(int m)
{
	assert( m > 0 );
	multiplicity_ = m;
}

std::vector<int> const &
Tetrahedron::get_corner_indices() const
{
	return cornerIndicesData_;
}

void
Tetrahedron::compute_corner_vectors(
		std::vector<double> & p0,
		std::vector<double> & v123 ) const
{
	std::vector<double> p0123(3*4);
	this->compute_corner_points(p0123);

	v123.resize(9);
	// construct the vectors v1 = p1-p0; v2 = p2-p0 ...
	for ( int iv = 1 ; iv < 4 ; ++iv)
		for ( int xi = 0 ; xi < 3 ; ++xi)
			v123[3*(iv-1)+xi] = p0123[3*iv+xi] - p0123[xi];
}
void
Tetrahedron::compute_corner_points(
		std::vector<double> & p0123 ) const
{
	assert(cornerIndicesExtendedReducible_.size() == 4);
	p0123.resize(12);
	for ( int iv = 0 ; iv < 4 ; ++iv)
	{
		auto red = cornerIndicesExtendedReducible_[iv];
		std::vector<double> p = extendedGrid_->get_vector_direct(red);
		for ( int xi = 0 ; xi < 3 ; ++xi)
			p0123[iv*3+xi] = p[xi];
	}
}

void
Tetrahedron::check_vectors_inside(
		std::vector<double> const & v,
		std::vector<bool> & inside,
		std::vector<double> & barycentricCoordinates) const
{
	assert(v.size()%3 == 0);
	int const numVects = v.size()/3;

	barycentricCoordinates.resize(4*numVects);
	inside.resize(numVects);

	if ( numVects == 0 )
		return;
	std::vector<double> p0123;
	this->compute_corner_points(p0123);

	double x1 = p0123[0];
	double y1 = p0123[1];
	double z1 = p0123[2];

	double x2 = p0123[3];
	double y2 = p0123[4];
	double z2 = p0123[5];

	double x3 = p0123[6];
	double y3 = p0123[7];
	double z3 = p0123[8];

	double x4 = p0123[9];
	double y4 = p0123[10];
	double z4 = p0123[11];

	// formula: http://steve.hollasch.net/cgindex/geometry/ptintet.html
	double d0 = -(x3*y2*z1) + x4*y2*z1 + x2*y3*z1 - x4*y3*z1 - x2*y4*z1 + x3*y4*z1 +
			   x3*y1*z2 - x4*y1*z2 - x1*y3*z2 + x4*y3*z2 + x1*y4*z2 - x3*y4*z2 -
			   x2*y1*z3 + x4*y1*z3 + x1*y2*z3 - x4*y2*z3 - x1*y4*z3 + x2*y4*z3 +
			   x2*y1*z4 - x3*y1*z4 - x1*y2*z4 + x3*y2*z4 + x1*y3*z4 - x2*y3*z4;

	for ( int ip = 0 ; ip < numVects ; ++ip )
	{
		double x = v[ip*3+0];
		double y = v[ip*3+1];
		double z = v[ip*3+2];

		double d1 = -(x3*y2*z) + x4*y2*z + x2*y3*z - x4*y3*z - x2*y4*z + x3*y4*z + x3*y*z2 -
				   x4*y*z2 - x*y3*z2 + x4*y3*z2 + x*y4*z2 - x3*y4*z2 - x2*y*z3 + x4*y*z3 +
				   x*y2*z3 - x4*y2*z3 - x*y4*z3 + x2*y4*z3 + x2*y*z4 - x3*y*z4 - x*y2*z4 +
				   x3*y2*z4 + x*y3*z4 - x2*y3*z4;
		double d2 = x3*y1*z - x4*y1*z - x1*y3*z + x4*y3*z + x1*y4*z - x3*y4*z - x3*y*z1 +
				   x4*y*z1 + x*y3*z1 - x4*y3*z1 - x*y4*z1 + x3*y4*z1 + x1*y*z3 - x4*y*z3 -
				   x*y1*z3 + x4*y1*z3 + x*y4*z3 - x1*y4*z3 - x1*y*z4 + x3*y*z4 + x*y1*z4 -
				   x3*y1*z4 - x*y3*z4 + x1*y3*z4;
		double d3 = -(x2*y1*z) + x4*y1*z + x1*y2*z - x4*y2*z - x1*y4*z + x2*y4*z + x2*y*z1 -
				   x4*y*z1 - x*y2*z1 + x4*y2*z1 + x*y4*z1 - x2*y4*z1 - x1*y*z2 + x4*y*z2 +
				   x*y1*z2 - x4*y1*z2 - x*y4*z2 + x1*y4*z2 + x1*y*z4 - x2*y*z4 - x*y1*z4 +
				   x2*y1*z4 + x*y2*z4 - x1*y2*z4;
		double d4 = x2*y1*z - x3*y1*z - x1*y2*z + x3*y2*z + x1*y3*z - x2*y3*z - x2*y*z1 +
				   x3*y*z1 + x*y2*z1 - x3*y2*z1 - x*y3*z1 + x2*y3*z1 + x1*y*z2 - x3*y*z2 -
				   x*y1*z2 + x3*y1*z2 + x*y3*z2 - x1*y3*z2 - x1*y*z3 + x2*y*z3 + x*y1*z3 -
				   x2*y1*z3 - x*y2*z3 + x1*y2*z3;

		assert((d0 - (d1+d2+d3+d4)) < 1e-8);

		// definition : point on the boundary is inside.
		// We allow a numerical grace zone what we consider "on" the border.
		// the value below is effectively a threshold on the volume fraction
		double signd0 = d0/std::abs(d0);
		const double borderThr = 1e-6;

		barycentricCoordinates[ip*4+0] = d1/d0;
		barycentricCoordinates[ip*4+1] = d2/d0;
		barycentricCoordinates[ip*4+2] = d3/d0;
		barycentricCoordinates[ip*4+3] = d4/d0;
		inside[ip] = std::all_of(&barycentricCoordinates[ip*4+0], &barycentricCoordinates[ip*4+0]+4,
				[&borderThr](double b){return (b > -borderThr) and (b < 1.0+borderThr); });
	}
}

bool operator< (Tetrahedron const & t1, Tetrahedron const & t2)
{
	assert(t1.cornerIndicesExtended_.size() == 4);
	assert(t2.cornerIndicesExtended_.size() == 4);
	for ( int i = 0 ; i < 4 ; ++i)
		if ( t1.cornerIndicesExtended_[i] != t2.cornerIndicesExtended_[i] )
			return t1.cornerIndicesExtended_[i] < t2.cornerIndicesExtended_[i];
	return false;
}

} /* namespace LatticeStructure */
} /* namespace elephon */
