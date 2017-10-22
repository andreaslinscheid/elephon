/*	This file TetrahedraGrid.cpp is part of elephon.
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
 *  Created on: Oct 8, 2017
 *      Author: A. Linscheid
 */

#include "LatticeStructure/TetrahedraGrid.h"
#include <algorithm>

namespace elephon
{
namespace LatticeStructure
{

void
TetrahedraGrid::initialize(std::shared_ptr<const RegularSymmetricGrid> grid)
{
	grid_ = grid;
	tetras_.clear();

	LatticeStructure::ExtendedSymmetricGrid extendedGrid;
	extendedGrid.initialize(
			{grid_->get_grid_dim()[0]+1, grid_->get_grid_dim()[1]+1, grid_->get_grid_dim()[2]+1},
			grid_->get_grid_prec(),
			grid_->get_grid_shift(),
			grid_->get_symmetry(),
			grid_->get_lattice());
	auto extendedGrid_ptr = std::make_shared<ExtendedSymmetricGrid>(std::move(extendedGrid));

	std::vector< RegularSymmetricGrid::GridCube > cubes;
	grid_->get_grid_cubes(cubes);

	auto norm2_indices = [&grid] (int index1, std::vector<double> const & d)
	{
		auto v1 = grid->get_vector_direct( index1 );
		for ( auto &xi : v1 )
			xi -= std::floor(xi);
		double diag2 = std::pow(v1[0]-d[0],2)+std::pow(v1[1]-d[1],2)+std::pow(v1[2]-d[2],2);
		return diag2;
	};

	// this method relies on the ordering chosen by
	// RegularBareGrid::compute_reducible_cube_indices_surrounding_nongrid_point
	RegularSymmetricGrid::GridCube const & c0 = cubes.at(0);
	std::vector<double> d_0_6{	1.0/double(grid_->get_grid_dim()[0]),
								1.0/double(grid_->get_grid_dim()[1]),
								1.0/double(grid_->get_grid_dim()[2])};
	std::vector<double> d_1_7{ -1.0/double(grid_->get_grid_dim()[0]),
								1.0/double(grid_->get_grid_dim()[1]),
								1.0/double(grid_->get_grid_dim()[2])};
	bool mainDiagon_0_6 = norm2_indices(c0.cornerIndices_[0], d_0_6)
						 < norm2_indices(c0.cornerIndices_[1], d_1_7);

	std::map<Tetrahedra,std::vector<std::int32_t>> tetraGrid;
	reducibleTetras_.clear();
	reducibleTetras_.reserve( this->get_n_reducible_tetra() );
	for ( auto const & c : cubes )
	{
		this->split_cube_insert_tetra(c, extendedGrid_ptr, mainDiagon_0_6, reducibleTetras_, tetraGrid);
	}

	tetras_.reserve( tetraGrid.size() );
	irreducibleToReducible_.reserve( tetraGrid.size() );
	reducibleToIrreducible_.assign( reducibleTetras_.size() , -1);
	for ( auto & p : tetraGrid )
	{
		for ( auto i : p.second)
		{
			assert( reducibleToIrreducible_[i] < 0 );
			reducibleToIrreducible_[i] = tetras_.size();
		}
		tetras_.push_back(p.first);
		tetras_.rbegin()->set_multiplicity(p.second.size());
		irreducibleToReducible_.push_back(std::move(p.second));
	}

	assert( *std::min_element(reducibleToIrreducible_.begin(), reducibleToIrreducible_.end()) == 0 );
}

void
TetrahedraGrid::split_cube_insert_tetra(
		RegularSymmetricGrid::GridCube const & cube,
		std::shared_ptr<const ExtendedSymmetricGrid> extendedGrid,
		bool diagonal1,
		std::vector<Tetrahedra> & reducibleTetra,
		std::map<Tetrahedra,std::vector<std::int32_t>> & tetraSet)
{
	auto xyz = grid_->get_reducible_to_xyz( cube.cornerIndices_[0] );
	std::vector<int> extendendCubeMap{
	 	extendedGrid->get_xyz_to_reducible(xyz),     // 1 = min, min, min
		extendedGrid->get_xyz_to_reducible({xyz[0]+1, xyz[1],   xyz[2]}),     // 2 = max, min, min
		extendedGrid->get_xyz_to_reducible({xyz[0]+1, xyz[1]+1, xyz[2]}),     // 3 = max, max, min
		extendedGrid->get_xyz_to_reducible({xyz[0],   xyz[1]+1, xyz[2]}),     // 4 = min, max, min
		extendedGrid->get_xyz_to_reducible({xyz[0],   xyz[1],   xyz[2]+1}),   // 5 = min, min, max
		extendedGrid->get_xyz_to_reducible({xyz[0]+1, xyz[1],   xyz[2]+1}),   // 6 = max, min, max
		extendedGrid->get_xyz_to_reducible({xyz[0]+1, xyz[1]+1, xyz[2]+1}),   // 7 = max, max, max
		extendedGrid->get_xyz_to_reducible({xyz[0],   xyz[1]+1, xyz[2]+1}),   // 8 = min, max, max
	};

	// this method relies on the ordering chosen by
	// RegularBareGrid::compute_reducible_cube_indices_surrounding_nongrid_point
	std::vector<std::vector<int>> tetraCornerReducibleIndexList(6);

	if ( diagonal1 )
	{
		tetraCornerReducibleIndexList[0] = std::vector<int>{0, 6, 4, 5};
		tetraCornerReducibleIndexList[1] = std::vector<int>{0, 6, 4, 7};
		tetraCornerReducibleIndexList[2] = std::vector<int>{0, 6, 5, 1};
		tetraCornerReducibleIndexList[3] = std::vector<int>{0, 6, 1, 2};
		tetraCornerReducibleIndexList[4] = std::vector<int>{0, 6, 2, 3};
		tetraCornerReducibleIndexList[5] = std::vector<int>{0, 6, 3, 7};
	}
	else
	{
		tetraCornerReducibleIndexList[0] = std::vector<int>{0, 1, 3, 5};
		tetraCornerReducibleIndexList[1] = std::vector<int>{1, 3, 2, 5};
		tetraCornerReducibleIndexList[2] = std::vector<int>{3, 5, 2, 6};
		tetraCornerReducibleIndexList[3] = std::vector<int>{3, 5, 7, 6};
		tetraCornerReducibleIndexList[4] = std::vector<int>{4, 5, 7, 3};
		tetraCornerReducibleIndexList[5] = std::vector<int>{0, 4, 5, 3};
	}

	for ( auto & t : tetraCornerReducibleIndexList)
	{
		int reducibleIndex = reducibleTetra.size();
		// here we convert indices within the cube to irregular grid indices in the extended
		// grid for tetrahedra identification.
		auto extendedIndices = t;
		for ( auto & index : extendedIndices)
			index = extendedGrid->get_maps_red_to_irreducible()[ extendendCubeMap[index] ];
		auto extendedIndicesReducible = t;
		for ( auto & index : extendedIndicesReducible)
			index = extendendCubeMap[index];

		auto reducibleDataIndices = t;
		for ( auto & index : reducibleDataIndices)
			index = cube.cornerIndices_[index];

		// here we convert indices within the cube to irregular grid indices for data lookup.
		for ( auto & index : t)
			index = grid_->get_maps_red_to_irreducible()[ cube.cornerIndices_[index] ];

		Tetrahedra reducibleTetrahedron(
				std::move(reducibleDataIndices),
				extendedIndicesReducible,
				extendedIndicesReducible,
				extendedGrid);
		reducibleTetra.push_back(std::move(reducibleTetrahedron));

		Tetrahedra tetra(
				std::move(t),
				std::move(extendedIndices),
				std::move(extendedIndicesReducible),
				extendedGrid);
		auto ret = tetraSet.insert( std::move(std::make_pair(tetra, std::vector<std::int32_t>())) );
		ret.first->second.push_back( reducibleIndex );
	}
}

int
TetrahedraGrid::get_n_tetra() const
{
	return tetras_.size();
}

int
TetrahedraGrid::get_n_reducible_tetra() const
{
	return 6*grid_->get_grid_dim()[0]*grid_->get_grid_dim()[1]*grid_->get_grid_dim()[2];
}

std::vector<TetrahedraGrid::Tetrahedra> const
TetrahedraGrid::get_tetra_list() const
{
	return tetras_;
}

std::vector<TetrahedraGrid::Tetrahedra> const
TetrahedraGrid::get_reducible_tetra_list() const
{
	return reducibleTetras_;
}

std::shared_ptr<const RegularSymmetricGrid>
TetrahedraGrid::get_grid() const
{
	return grid_;
}

int
TetrahedraGrid::get_reducible_to_irreducible(int ired) const
{
	assert((ired>=0)&&(ired<reducibleToIrreducible_.size()));
	return reducibleToIrreducible_[ired];
}

std::vector<int> const &
TetrahedraGrid::get_irreducible_to_reducible(int iirred) const
{
	assert((iirred>=0)&&(iirred<irreducibleToReducible_.size()));
	return irreducibleToReducible_[iirred];
}


namespace detail
{
Tetrahedra::Tetrahedra(
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
Tetrahedra::get_multiplicity() const
{
	return multiplicity_;
}

void
Tetrahedra::set_multiplicity(int m)
{
	assert( m > 0 );
	multiplicity_ = m;
}

std::vector<int> const &
Tetrahedra::get_corner_indices() const
{
	return cornerIndicesData_;
}

void
Tetrahedra::compute_corner_vectors(
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
Tetrahedra::compute_corner_points(
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

bool operator< (Tetrahedra const & t1, Tetrahedra const & t2)
{
	assert(t1.cornerIndicesExtended_.size() == 4);
	assert(t2.cornerIndicesExtended_.size() == 4);
	for ( int i = 0 ; i < 4 ; ++i)
		if ( t1.cornerIndicesExtended_[i] != t2.cornerIndicesExtended_[i] )
			return t1.cornerIndicesExtended_[i] < t2.cornerIndicesExtended_[i];
	return false;
}
} /* namespace detail */

} /* namespace LatticeStructure */
} /* namespace elephon */
