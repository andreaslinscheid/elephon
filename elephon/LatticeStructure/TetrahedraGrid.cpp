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
#include <map>
#include <set>
#include <iostream>

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
	extendedGrid_ = std::make_shared<ExtendedSymmetricGrid>(std::move(extendedGrid));

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
	mainDiagon_0_6_ = norm2_indices(c0.cornerIndices_[0], d_0_6)
						 < norm2_indices(c0.cornerIndices_[1], d_1_7);

	std::map<Tetrahedron,std::vector<std::int32_t>> tetraGrid;
	reducibleTetras_.clear();
	reducibleTetras_.reserve( this->get_n_reducible_tetra() );
	for ( auto const & c : cubes )
	{
		this->split_cube_insert_tetra(c, reducibleTetras_, tetraGrid);
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
		std::vector<Tetrahedron> & reducibleTetra,
		std::map<Tetrahedron,std::vector<std::int32_t>> & tetraSet) const
{
	auto xyz = grid_->get_reducible_to_xyz( cube.cornerIndices_[0] );
	std::vector<int> extendendCubeMap{
		extendedGrid_->get_xyz_to_reducible(xyz),     // 1 = min, min, min
		extendedGrid_->get_xyz_to_reducible({xyz[0]+1, xyz[1],   xyz[2]}),     // 2 = max, min, min
		extendedGrid_->get_xyz_to_reducible({xyz[0]+1, xyz[1]+1, xyz[2]}),     // 3 = max, max, min
		extendedGrid_->get_xyz_to_reducible({xyz[0],   xyz[1]+1, xyz[2]}),     // 4 = min, max, min
		extendedGrid_->get_xyz_to_reducible({xyz[0],   xyz[1],   xyz[2]+1}),   // 5 = min, min, max
		extendedGrid_->get_xyz_to_reducible({xyz[0]+1, xyz[1],   xyz[2]+1}),   // 6 = max, min, max
		extendedGrid_->get_xyz_to_reducible({xyz[0]+1, xyz[1]+1, xyz[2]+1}),   // 7 = max, max, max
		extendedGrid_->get_xyz_to_reducible({xyz[0],   xyz[1]+1, xyz[2]+1}),   // 8 = min, max, max
	};

	// this method relies on the ordering chosen by
	// RegularBareGrid::compute_reducible_cube_indices_surrounding_nongrid_point
	std::vector<std::vector<int>> tetraCornerReducibleIndexList(6);

	if ( mainDiagon_0_6_ )
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
			index = extendedGrid_->get_maps_red_to_irreducible()[ extendendCubeMap[index] ];
		auto extendedIndicesReducible = t;
		for ( auto & index : extendedIndicesReducible)
			index = extendendCubeMap[index];

		auto reducibleDataIndices = t;
		for ( auto & index : reducibleDataIndices)
			index = cube.cornerIndices_[index];

		// here we convert indices within the cube to irregular grid indices for data lookup.
		for ( auto & index : t)
			index = grid_->get_maps_red_to_irreducible()[ cube.cornerIndices_[index] ];

		Tetrahedron reducibleTetrahedron(
				std::move(reducibleDataIndices),
				extendedIndicesReducible,
				extendedIndicesReducible,
				extendedGrid_);
		reducibleTetra.push_back(std::move(reducibleTetrahedron));

		Tetrahedron tetra(
				std::move(t),
				std::move(extendedIndices),
				std::move(extendedIndicesReducible),
				extendedGrid_);
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

std::vector<Tetrahedron> const
TetrahedraGrid::get_tetra_list() const
{
	return tetras_;
}

std::vector<Tetrahedron> const
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

void
TetrahedraGrid::compute_grid_tetrahedra_surrounding_nongrid_points(
		std::vector<double> const & nonGridPoints,
		std::map<Tetrahedron,std::vector<int>> & tetras) const
{
	std::vector<int> nonGridPtToCubeMap;
	std::vector<RegularBareGrid::GridCube> cubes;
	grid_->compute_grid_cubes_surrounding_nongrid_points(
			nonGridPoints,
			nonGridPtToCubeMap,
			cubes);

	tetras.clear();
	std::vector<bool> insideTetra;
	std::vector<double> vectors, barycen;
	for ( auto const & cube : cubes )
	{
		// build all 6 tetrahedra
		std::vector<Tetrahedron> reducibleTetraThisCube;
		std::map<Tetrahedron,std::vector<std::int32_t>> dummy;
		this->split_cube_insert_tetra(cube, reducibleTetraThisCube, dummy);

		// see for each contained grid point in which one it is
		// then add the respective tetrahedron to tetras, appending
		// the grid vector index to the vector pointed to.
		// Since tetrahedra evaluate a grid point on the border as true
		// we need to keep track of the points we distribute so that no
		// point appears twice
		std::set<int> availableGridIndices(cube.containedIrregularPts_.begin(), cube.containedIrregularPts_.end());
		int numAvail = cube.containedIrregularPts_.size();
		for (auto const & t : reducibleTetraThisCube)
		{
			if (availableGridIndices.size() == 0)
				break;

			// copy the vectors in question and map them to the zone [0,1[
			std::vector<int> availableGridIndicesVector(availableGridIndices.begin(), availableGridIndices.end());
			vectors.resize(availableGridIndices.size()*3);
			for (int pi  = 0 ; pi < availableGridIndicesVector.size(); ++pi)
				for (int i = 0 ; i < 3 ; ++i)
					vectors[pi*3+i] = nonGridPoints[availableGridIndicesVector[pi]*3+i]
									-std::floor(nonGridPoints[availableGridIndicesVector[pi]*3+i]);

			std::vector<int> containedNonGridPoints;
			t.check_vectors_inside(vectors, insideTetra, barycen);
			assert(insideTetra.size() == availableGridIndicesVector.size());
			for ( int ip = 0 ; ip < availableGridIndicesVector.size() ; ++ip)
			{
				if ( insideTetra[ip] )
				{
					containedNonGridPoints.push_back(availableGridIndicesVector[ip]);
					// remove the indices from the set so that they are not present in the next tetrahedron search
					// this way, a given point will be only part of one tetrahedron, even if that point is on the border.
					availableGridIndices.erase(availableGridIndicesVector[ip]);
					numAvail--;
				}
			}
			if ( not containedNonGridPoints.empty() )
			{
				auto it = tetras.insert(std::make_pair(t, std::vector<int>()));
				it.first->second = std::move(containedNonGridPoints);
			}
		}

		if ( numAvail != 0 )
		{
			std::cout << "The follwing vectors were not assigned ";
			for (auto i : availableGridIndices)
			{
				std::vector<double> bla{ nonGridPoints[i*3+0],  nonGridPoints[i*3+1],  nonGridPoints[i*3+2]};
				std::cout << '\n' << bla[0]<<'\t'<< bla[1]<<'\t'<< bla[2];
				reducibleTetraThisCube[0].check_vectors_inside(bla, insideTetra, barycen);
				std::cout << '\n'<< insideTetra[0];
				std::cout << '\n' << barycen[0*4+0]<<'\t'<< barycen[0*4+1]<<'\t'<< barycen[0*4+2];
				std::cout << '\n' << barycen[1*4+0]<<'\t'<< barycen[1*4+1]<<'\t'<< barycen[1*4+2];
				std::cout << '\n' << barycen[2*4+0]<<'\t'<< barycen[2*4+1]<<'\t'<< barycen[2*4+2];
				std::cout << '\n' << barycen[3*4+0]<<'\t'<< barycen[3*4+1]<<'\t'<< barycen[3*4+2];


				reducibleTetraThisCube[1].check_vectors_inside(bla, insideTetra, barycen);
				std::cout << "\n\n"<< insideTetra[0];
				std::cout << '\n' << barycen[0*4+0]<<'\t'<< barycen[0*4+1]<<'\t'<< barycen[0*4+2];
				std::cout << '\n' << barycen[1*4+0]<<'\t'<< barycen[1*4+1]<<'\t'<< barycen[1*4+2];
				std::cout << '\n' << barycen[2*4+0]<<'\t'<< barycen[2*4+1]<<'\t'<< barycen[2*4+2];
				std::cout << '\n' << barycen[3*4+0]<<'\t'<< barycen[3*4+1]<<'\t'<< barycen[3*4+2];

				reducibleTetraThisCube[2].check_vectors_inside(bla, insideTetra, barycen);
				std::cout << "\n\n"<< insideTetra[0];
				std::cout << '\n' << barycen[0*4+0]<<'\t'<< barycen[0*4+1]<<'\t'<< barycen[0*4+2];
				std::cout << '\n' << barycen[1*4+0]<<'\t'<< barycen[1*4+1]<<'\t'<< barycen[1*4+2];
				std::cout << '\n' << barycen[2*4+0]<<'\t'<< barycen[2*4+1]<<'\t'<< barycen[2*4+2];
				std::cout << '\n' << barycen[3*4+0]<<'\t'<< barycen[3*4+1]<<'\t'<< barycen[3*4+2];

				reducibleTetraThisCube[3].check_vectors_inside(bla, insideTetra, barycen);
				std::cout << "\n\n"<< insideTetra[0];
				std::cout << '\n' << barycen[0*4+0]<<'\t'<< barycen[0*4+1]<<'\t'<< barycen[0*4+2];
				std::cout << '\n' << barycen[1*4+0]<<'\t'<< barycen[1*4+1]<<'\t'<< barycen[1*4+2];
				std::cout << '\n' << barycen[2*4+0]<<'\t'<< barycen[2*4+1]<<'\t'<< barycen[2*4+2];
				std::cout << '\n' << barycen[3*4+0]<<'\t'<< barycen[3*4+1]<<'\t'<< barycen[3*4+2];

				reducibleTetraThisCube[4].check_vectors_inside(bla, insideTetra, barycen);
				std::cout << "\n\n"<< insideTetra[0];
				std::cout << '\n' << barycen[0*4+0]<<'\t'<< barycen[0*4+1]<<'\t'<< barycen[0*4+2];
				std::cout << '\n' << barycen[1*4+0]<<'\t'<< barycen[1*4+1]<<'\t'<< barycen[1*4+2];
				std::cout << '\n' << barycen[2*4+0]<<'\t'<< barycen[2*4+1]<<'\t'<< barycen[2*4+2];
				std::cout << '\n' << barycen[3*4+0]<<'\t'<< barycen[3*4+1]<<'\t'<< barycen[3*4+2];

				reducibleTetraThisCube[5].check_vectors_inside(bla, insideTetra, barycen);
				std::cout << "\n\n"<< insideTetra[0];
				std::cout << '\n' << barycen[0*4+0]<<'\t'<< barycen[0*4+1]<<'\t'<< barycen[0*4+2];
				std::cout << '\n' << barycen[1*4+0]<<'\t'<< barycen[1*4+1]<<'\t'<< barycen[1*4+2];
				std::cout << '\n' << barycen[2*4+0]<<'\t'<< barycen[2*4+1]<<'\t'<< barycen[2*4+2];
				std::cout << '\n' << barycen[3*4+0]<<'\t'<< barycen[3*4+1]<<'\t'<< barycen[3*4+2];
			}
			std::cout << std::endl;
			throw std::logic_error("Did not distribute all points within a cube"
					" to the 6 tetrahedra that completely fill the cube?");
		}
	}
}

} /* namespace LatticeStructure */
} /* namespace elephon */
