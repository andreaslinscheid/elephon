/*	This file TetrahedraIsosurface.cpp is part of elephon.
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
 *  Created on: Oct 9, 2017
 *      Author: A. Linscheid
 */

#include "ElectronicStructure/TetrahedraIsosurface.h"
#include "Algorithms/helperfunctions.hpp"
#include <vtkPoints.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkTriangle.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>

namespace elephon
{
namespace ElectronicStructure
{

namespace detail
{
struct PrelimEdge
{
	PrelimEdge(int p1, int p2) {
		p1_ = p1 < p2 ? p1 : p2;
		p2_ = p1 < p2 ? p2 : p1;
	};

	int p1_, p2_;

	mutable std::vector<int> triangleIndIices_;

	mutable std::set<int> tetraIndIices_;

	bool operator< (PrelimEdge const & e) const {
		if ( (not (e.p1_ < p1_)) and (not (p1_ < e.p1_ )) ) // p1 == e.p1_
			return p2_ < e.p2_;
		return p1_ < e.p1_;
	}
} ;

struct PrelimTriangle{
	PrelimTriangle(NonGridPoint p1, NonGridPoint p2, NonGridPoint p3, double weight, int multiplicity)
		: p1_(p1), p2_(p2), p3_(p3), weight_(weight), multiplicity_(multiplicity){};
	NonGridPoint p1_, p2_, p3_;
	double weight_;
	int multiplicity_ = 0;
};

struct PrelimTetra{
	PrelimTetra(NonGridPoint p1, NonGridPoint p2, NonGridPoint p3, NonGridPoint p4, double weight, int multiplicity)
		: p1_(p1), p2_(p2), p3_(p3),  p4_(p4), weight_(weight), multiplicity_(multiplicity){};
	NonGridPoint p1_, p2_, p3_, p4_;
	double weight_;
	int multiplicity_ = 0;
};

struct kinfo{
	double weight = 0;
	double reducible_weight = 0;
	std::vector<int> triangleIndex;
	int pointIndex = -1;
	int multiplicity = 0;
};
} /* namespace detail */

void
TetrahedraIsosurface::initialize(
		std::shared_ptr<const LatticeStructure::TetrahedraGrid> tetra,
		std::shared_ptr<const LatticeStructure::DataRegularGrid<double>> bands,
		std::vector<double> const & isoEnergies)
{
	assert( tetra );
	tetra_ = tetra;
	assert( bands );

	const double eps = bands->get_grid().get_symmetry().get_symmetry_prec();
	numBands_ = bands->get_nData_gpt();
	numIsoE_ = isoEnergies.size();

	std::vector<std::vector<detail::PrelimTriangle>> prlimTrgls(numBands_*numIsoE_);
	std::vector<std::vector<detail::PrelimTetra>> prlimTetras(numBands_*numIsoE_);

	std::vector<double> kIsos;
	for ( auto const & t : tetra->get_tetra_list() )
	{
		double tw = double(t.get_multiplicity()) / double(tetra->get_n_reducible_tetra());
		for ( int ib = 0 ; ib < numBands_ ; ++ib )
		{
			std::multimap<double,int> ecorner;
			for (int ie = 0 ; ie < 4 ; ++ie )
				ecorner.insert(std::make_pair(bands->read(t.get_corner_indices()[ie], ib),ie));

			assert(ecorner.size() == 4);
			auto it = ecorner.begin();
			double e1 = it->first;
			++it;
			double e2 = it->first;
			++it;
			double e3 = it->first;
			++it;
			double e4 = it->first;

			for ( int iw = 0 ; iw < numIsoE_ ; ++iw)
			{
				auto e = isoEnergies[iw];
				// formula is the energy derivative of 10.1103/PhysRevB.49.16223 Appendix A
				if ( (e < e1) or ( e > e4 ) )
					continue;

				double kweight = 0;

				if ( (e >= e1) and (e < e2) )
					kweight = tw*(3.0*std::pow(e-e1,2))/(e2-e1)/(e3-e1)/(e4-e1);

				if ( (e >= e2) and (e < e3) )
					kweight = tw / (e3-e1)/(e4-e1)
								*(3*(e2-e1) + 3*2*(e-e2) - ((e3-e1)+(e4-e2))/(e3-e2)/(e4-e2)*3*std::pow(e-e2,2) );

				if ( (e >= e3) and (e <= e4) )
					kweight = tw*3*std::pow(e4-e,2)/(e4-e1)/(e4-e2)/(e4-e3);

				this->compute_iso_vectors(t, ecorner, e, kIsos);

				int nkIsoAlgo = kIsos.size()/3;
				assert((nkIsoAlgo == 4) or (nkIsoAlgo == 3));
				// first check if points are actually equivalent up to grid precision. this is automatically done
				// by using this set
				std::set<NonGridPoint> newPoints;
				for ( int inp = 0 ; inp < nkIsoAlgo ; ++inp )
					newPoints.insert(NonGridPoint(kIsos[inp*3+0], kIsos[inp*3+1], kIsos[inp*3+2], eps));

				// it can happen that we arrive at down to a single point. In cases where we end up with < 3 points
				// we create empty "triangles" with vertices that point to the same point.
				int nkIso = newPoints.size();
				assert(nkIso > 0);

				if (nkIso == 1 )
				{
					auto ngp = *newPoints.begin();
					prlimTrgls[iw*numBands_+ib].push_back(
							detail::PrelimTriangle(ngp, ngp, ngp, kweight, t.get_multiplicity()) );
				}

				if (nkIso == 2 )
				{
					auto ngp1 = *newPoints.begin();
					auto ngp2 = *(newPoints.begin()++);
					prlimTrgls[iw*numBands_+ib].push_back(
							detail::PrelimTriangle(ngp1, ngp2, ngp1, kweight, t.get_multiplicity()) );
				}

				// extract connectivity information. In case of one triangle this is easy => all k point
				// belong to the newly added single triangle.
				if (nkIso == 3 )
				{
					std::vector<NonGridPoint> pts;
					for ( int ikIso = 0 ; ikIso < nkIso; ++ikIso )
						pts.push_back(NonGridPoint(kIsos[ikIso*3+0], kIsos[ikIso*3+1], kIsos[ikIso*3+2], eps));
					prlimTrgls[iw*numBands_+ib].push_back(
							detail::PrelimTriangle(pts[0], pts[1], pts[2], kweight, t.get_multiplicity()) );
				}

				if (nkIso == 4 )
				{
					// decide how to split the 4 points into two triangles. This is
					// not an easy problem. We shall insert all 4 sides of the
					// tetrahedron and remove two later in a way such that we do
					// not create a whole.
					auto fetchPoint = [&kIsos] (int i){ return std::vector<double>(&kIsos[i*3], &kIsos[i*3]+3);};
					auto p1 = fetchPoint(0);
					auto p2 = fetchPoint(1);
					auto p3 = fetchPoint(2);
					auto p4 = fetchPoint(3);
					NonGridPoint ngp1(p1[0], p1[1], p1[2], eps);
					NonGridPoint ngp2(p2[0], p2[1], p2[2], eps);
					NonGridPoint ngp3(p3[0], p3[1], p3[2], eps);
					NonGridPoint ngp4(p4[0], p4[1], p4[2], eps);
					prlimTetras[iw*numBands_+ib].push_back(
							detail::PrelimTetra(ngp1, ngp2, ngp3, ngp4, kweight, t.get_multiplicity()) );
				}
			}
		}
	}

	// now, we need to use the points stored in the triangles to connect the mesh.
	// this means we have to collapse equal points. A problem is that a case can
	// happen where p1 == p2 && p2 == p3 but p1 != p3 due to the finite cutoff.
	// in that case the order of insertion matters and we can end up with gaps in the mesh.
	kIso_.resize(numBands_*numIsoE_);
	kIsoWeights_.resize(numBands_*numIsoE_);
	kIsoReducibleWeights_.resize(numBands_*numIsoE_);
	triangles_.resize(numBands_*numIsoE_);
	std::vector<std::map<NonGridPoint, detail::kinfo>> isoKandWeights(numBands_*numIsoE_);
	for ( int iw = 0 ; iw < numIsoE_ ; ++iw)
		for ( int ib = 0 ; ib < numBands_ ; ++ib )
		{
			std::set<detail::PrelimEdge> edges;
			int numTriangles = prlimTrgls[iw*numBands_+ib].size() + prlimTetras[iw*numBands_+ib].size()*2;
			auto & bandIsoSurface = isoKandWeights[iw*numBands_+ib];
			std::vector<int> triangleToPointIndexMap;
			triangleToPointIndexMap.reserve( numTriangles*3 );

			auto insertpt = [&bandIsoSurface] (NonGridPoint const & p, int itriangle, double weight, int objMultiplicity) {
				assert((p.k_[0] == p.k_[0]) and (p.k_[1] == p.k_[1]) and (p.k_[2] == p.k_[2]));
				auto ret = bandIsoSurface.insert( std::make_pair(p, detail::kinfo()) );
				ret.first->second.weight += weight;
				ret.first->second.reducible_weight += weight/objMultiplicity;
				if ( itriangle >= 0 )
					ret.first->second.triangleIndex.push_back(itriangle);
				if ( ret.first->second.pointIndex < 0)
					ret.first->second.pointIndex = static_cast<int>(bandIsoSurface.size())-1;
				return ret.first->second.pointIndex;
			};

			auto insert_triangle = [&] (
					NonGridPoint const & p1, NonGridPoint const & p2, NonGridPoint const & p3, double weight, int objMultiplicity)
				{
				int index = triangleToPointIndexMap.size()/3;
				int ip1 = insertpt(p1, index , weight, objMultiplicity);
				int ip2 = insertpt(p2, index , weight, objMultiplicity);
				int ip3 = insertpt(p3, index , weight, objMultiplicity);
				detail::PrelimEdge edge1(ip1, ip2);
				detail::PrelimEdge edge2(ip1, ip3);
				detail::PrelimEdge edge3(ip3, ip2);
				auto ret = edges.insert( edge1 );
				ret.first->triangleIndIices_.push_back(index);
				ret = edges.insert( edge2 );
				ret.first->triangleIndIices_.push_back(index);
				ret = edges.insert( edge3 );
				ret.first->triangleIndIices_.push_back(index);
				triangleToPointIndexMap.insert(triangleToPointIndexMap.end(), {ip1, ip2, ip3});
				};


			// first we add the simple triangles to the mesh.
			for (int it = 0 ; it < prlimTrgls[iw*numBands_+ib].size(); ++it)
			{
				auto const & tprlm = prlimTrgls[iw*numBands_+ib][it];
				insert_triangle(tprlm.p1_, tprlm.p2_, tprlm.p3_, tprlm.weight_/3.0, tprlm.multiplicity_);
			}

			// now we add the tetrahedra
			for (int it = 0 ; it < prlimTetras[iw*numBands_+ib].size(); ++it)
			{
				auto const & tprlm = prlimTetras[iw*numBands_+ib][it];
				int p1 = insertpt(tprlm.p1_, -1 , tprlm.weight_/4.0, tprlm.multiplicity_);
				int p2 = insertpt(tprlm.p2_, -1 , tprlm.weight_/4.0, tprlm.multiplicity_);
				int p3 = insertpt(tprlm.p3_, -1 , tprlm.weight_/4.0, tprlm.multiplicity_);
				int p4 = insertpt(tprlm.p4_, -1 , tprlm.weight_/4.0, tprlm.multiplicity_);
				auto ret = edges.insert( detail::PrelimEdge(p1, p2) );
				ret.first->tetraIndIices_.insert(it);
				ret = edges.insert( detail::PrelimEdge(p1, p3) );
				ret.first->tetraIndIices_.insert(it);
				ret = edges.insert( detail::PrelimEdge(p3, p2) );
				ret.first->tetraIndIices_.insert(it);
				ret = edges.insert( detail::PrelimEdge(p1, p4) );
				ret.first->tetraIndIices_.insert(it);
				ret = edges.insert( detail::PrelimEdge(p4, p3) );
				ret.first->tetraIndIices_.insert(it);
				ret = edges.insert( detail::PrelimEdge(p4, p2) );
				ret.first->tetraIndIices_.insert(it);
			}

			// now all points are correctly inserted. To get the triangulation bonus,
			// we have to discard the 2 of the 4 of each isotetrahedra. To not create
			// holes in the mesh, we do not discard edges that share any other two triangle
			// or tetrahedra. For rim edges, we discard one triangle of more than 2 sharing an edge
			auto get_edge = [&edges, &bandIsoSurface] (
					NonGridPoint const & p1,
					NonGridPoint const & p2,
					int tetraIndex) {
				auto p1_itr = bandIsoSurface.find(p1);
				auto p2_itr = bandIsoSurface.find(p2);
				assert( (p1_itr != bandIsoSurface.end()) && (p2_itr != bandIsoSurface.end()));
				auto it = edges.find(
						detail::PrelimEdge(p1_itr->second.pointIndex, p2_itr->second.pointIndex) );
				assert(it != edges.end());
				assert(it->tetraIndIices_.find(tetraIndex) != it->tetraIndIices_.end());
				it->tetraIndIices_.erase(tetraIndex);
				return it;
			};

			for (int it = 0 ; it < prlimTetras[iw*numBands_+ib].size(); ++it)
			{
				auto const & tprlm = prlimTetras[iw*numBands_+ib][it];
				auto e_it_1 = get_edge( tprlm.p1_, tprlm.p2_, it);
				auto e_it_2 = get_edge( tprlm.p1_, tprlm.p3_, it);
				auto e_it_3 = get_edge( tprlm.p3_, tprlm.p2_, it);
				auto e_it_4 = get_edge( tprlm.p4_, tprlm.p2_, it);
				auto e_it_5 = get_edge( tprlm.p1_, tprlm.p4_, it);
				auto e_it_6 = get_edge( tprlm.p3_, tprlm.p4_, it);

				auto get_neigb = [] (decltype(e_it_1) it )
					{
					return int(it->tetraIndIices_.size()+it->triangleIndIices_.size());
					};

				int numNeighbors = get_neigb(e_it_1) +  get_neigb(e_it_2) +
						get_neigb(e_it_3) +  get_neigb(e_it_4) +
						get_neigb(e_it_5) +  get_neigb(e_it_6);
				assert((numNeighbors >= 0) && (numNeighbors <= 4));

				if ( numNeighbors == 4 )
				{
					// an interior iso-tetrahedra has 4 neighbors sharing edges.
					// remove the first edge which does not have a neighbor
					if ( get_neigb(e_it_1) == 0 ) // p1 <-> p2 is cut
					{
						insert_triangle( tprlm.p1_, tprlm.p3_, tprlm.p4_, 0, tprlm.multiplicity_);
						insert_triangle( tprlm.p2_, tprlm.p3_, tprlm.p4_, 0, tprlm.multiplicity_);
						continue;
					}
					else if ( get_neigb(e_it_2) == 0 ) // p1 <-> p3 is cut
					{
						insert_triangle( tprlm.p1_, tprlm.p2_, tprlm.p4_, 0, tprlm.multiplicity_);
						insert_triangle( tprlm.p3_, tprlm.p2_, tprlm.p4_, 0, tprlm.multiplicity_);
						continue;
					}
					else if ( get_neigb(e_it_3) == 0 ) // p2 <-> p3 is cut
					{
						insert_triangle( tprlm.p2_, tprlm.p1_, tprlm.p4_, 0, tprlm.multiplicity_);
						insert_triangle( tprlm.p3_, tprlm.p1_, tprlm.p4_, 0, tprlm.multiplicity_);
						continue;
					}
					else if ( get_neigb(e_it_4) == 0 ) // p4 <-> p2 is cut
					{
						insert_triangle( tprlm.p2_, tprlm.p1_, tprlm.p3_, 0, tprlm.multiplicity_);
						insert_triangle( tprlm.p4_, tprlm.p1_, tprlm.p3_, 0, tprlm.multiplicity_);
						continue;
					}
					else if ( get_neigb(e_it_5) == 0 ) // p4 <-> p1 is cut
					{
						insert_triangle( tprlm.p1_, tprlm.p2_, tprlm.p3_, 0, tprlm.multiplicity_);
						insert_triangle( tprlm.p4_, tprlm.p2_, tprlm.p3_, 0, tprlm.multiplicity_);
						continue;
					}
					else if ( get_neigb(e_it_6) == 0 ) // p4 <-> p3 is cut
					{
						insert_triangle( tprlm.p3_, tprlm.p2_, tprlm.p1_, 0, tprlm.multiplicity_);
						insert_triangle( tprlm.p4_, tprlm.p2_, tprlm.p1_, 0, tprlm.multiplicity_);
						continue;
					}
					throw std::logic_error("Incorrect neighbor determination");
				}

				if ( numNeighbors == 3 )
				{
					// this is a rim iso-tetrahedron. A correct triangluation does create triangles
					// that share an edge within the iso-tetrahedron. In a first step, we insert a
					// valid first triangle with two edges shared by neighbors and then in the next
					// step, we only have one choice given that the number of neighbors must be <3.
					if (((get_neigb(e_it_1) == 1) and (get_neigb(e_it_2) == 1)) or
						((get_neigb(e_it_1) == 1) and (get_neigb(e_it_3) == 1)) or
						((get_neigb(e_it_2) == 1) and (get_neigb(e_it_3) == 1)) )
					{
						insert_triangle( tprlm.p1_, tprlm.p2_, tprlm.p3_,0, tprlm.multiplicity_);
					}
					else if (((get_neigb(e_it_1) == 1) and (get_neigb(e_it_4) == 1)) or
							 ((get_neigb(e_it_1) == 1) and (get_neigb(e_it_5) == 1)) or
							 ((get_neigb(e_it_4) == 1) and (get_neigb(e_it_5) == 1)) )
					{
						insert_triangle( tprlm.p1_, tprlm.p2_, tprlm.p4_, 0, tprlm.multiplicity_);
					}
					else if (((get_neigb(e_it_2) == 1) and (get_neigb(e_it_5) == 1)) or
							 ((get_neigb(e_it_2) == 1) and (get_neigb(e_it_6) == 1)) or
							 ((get_neigb(e_it_5) == 1) and (get_neigb(e_it_6) == 1)) )
					{
						insert_triangle( tprlm.p1_, tprlm.p3_, tprlm.p4_, 0, tprlm.multiplicity_);
					}
					else if (((get_neigb(e_it_4) == 1) and (get_neigb(e_it_3) == 1)) or
							 ((get_neigb(e_it_6) == 1) and (get_neigb(e_it_3) == 1)) or
							 ((get_neigb(e_it_6) == 1) and (get_neigb(e_it_4) == 1)) )
					{
						insert_triangle( tprlm.p3_, tprlm.p2_, tprlm.p4_, 0, tprlm.multiplicity_);
					}

					if ((get_neigb(e_it_1) != 2) and (get_neigb(e_it_2) != 2) and (get_neigb(e_it_3) != 2))
					{
						insert_triangle( tprlm.p1_, tprlm.p2_, tprlm.p3_, 0, tprlm.multiplicity_);
						continue;
					}
					else if ((get_neigb(e_it_1) != 2) and (get_neigb(e_it_4) != 2) and (get_neigb(e_it_5) != 2))
					{
						insert_triangle( tprlm.p1_, tprlm.p2_, tprlm.p4_, 0, tprlm.multiplicity_);
						continue;
					}
					else if ((get_neigb(e_it_2) != 2) and (get_neigb(e_it_5) != 2) and (get_neigb(e_it_6) != 2))
					{
						insert_triangle( tprlm.p1_, tprlm.p3_, tprlm.p4_, 0, tprlm.multiplicity_);
						continue;
					}
					else if ((get_neigb(e_it_3) != 2) and (get_neigb(e_it_4) != 2) and (get_neigb(e_it_6) != 2))
					{
						insert_triangle( tprlm.p3_, tprlm.p2_, tprlm.p4_, 0, tprlm.multiplicity_);
						continue;
					}
					throw std::logic_error("Incorrect neighbor determination neighbors = 3");
				}

				if ( numNeighbors == 2 )
				{
					// This can be an rim isotriangle on a corner if the edges with neighbors
					// share a common point.
					std::map<int, std::vector<decltype(e_it_1)> > pointIndices;
					auto check_edge = [&] (decltype(e_it_1) edgeIterator) {
						if (get_neigb(edgeIterator) == 1)
						{
							auto ret = pointIndices.insert(
									std::make_pair(edgeIterator->p1_, std::vector<decltype(e_it_1)>() ) );
							ret.first->second.push_back(edgeIterator);
							ret = pointIndices.insert(
									std::make_pair(edgeIterator->p2_, std::vector<decltype(e_it_1)>() ) );
							ret.first->second.push_back(edgeIterator);
						}
					};
					check_edge(e_it_1);
					check_edge(e_it_2);
					check_edge(e_it_3);
					check_edge(e_it_4);
					check_edge(e_it_5);
					check_edge(e_it_6);
					assert( pointIndices.size() > 2 );

					// If we face 4 points the iso tetrahedron connects two distinct pieces of
					// surface. This is difficult to resolve, since it cannot be easily determined
					// which orientation is "right". We simply add 2 triangles.
					if (pointIndices.size() == 4 )
					{
						insert_triangle( tprlm.p1_, tprlm.p3_, tprlm.p4_, 0, tprlm.multiplicity_);
						insert_triangle( tprlm.p2_, tprlm.p3_, tprlm.p4_, 0, tprlm.multiplicity_);
						continue;
					}
					else if (pointIndices.size() == 3 )
					{
						// We have two edges sharing a common point. We can obtain a valid triangulation
						// by choosing two triangles that share a common edge which is not one with
						// neighbors and have another edge with a neighbor each.
						// first we find the common point and the two edges with neighbors
						decltype(e_it_1) edge1 , edge2 ;
						auto it = pointIndices.begin();
						while( not (it->second.size() > 1) )
						{
							++it;
							if (it == pointIndices.end() )
								throw std::logic_error("Problem locating shared point of rim edge isotriangle");
						}
						edge1 = it->second[0];
						edge2 = it->second[1];

						// p1 shared point, p1->p4 shared edge
						if ( (edge1 == e_it_1) and (edge2 == e_it_2) )
						{
							insert_triangle( tprlm.p1_, tprlm.p4_, tprlm.p2_, 0, tprlm.multiplicity_);
							insert_triangle( tprlm.p1_, tprlm.p4_, tprlm.p3_, 0, tprlm.multiplicity_);
							continue;
						}

						// p2 shared point, p2->p4 shared edge
						if ( (edge1 == e_it_1) and (edge2 == e_it_3) )
						{
							insert_triangle( tprlm.p2_, tprlm.p4_, tprlm.p1_, 0, tprlm.multiplicity_);
							insert_triangle( tprlm.p2_, tprlm.p4_, tprlm.p3_, 0, tprlm.multiplicity_);
							continue;
						}

						// p2 shared point, p2->p3 shared edge
						if ( (edge1 == e_it_1) and (edge2 == e_it_4) )
						{
							insert_triangle( tprlm.p2_, tprlm.p3_, tprlm.p1_, 0, tprlm.multiplicity_);
							insert_triangle( tprlm.p2_, tprlm.p3_, tprlm.p4_, 0, tprlm.multiplicity_);
							continue;
						}

						// p1 shared point, p1->p3 shared edge
						if ( (edge1 == e_it_1) and (edge2 == e_it_5) )
						{
							insert_triangle( tprlm.p1_, tprlm.p3_, tprlm.p2_, 0, tprlm.multiplicity_);
							insert_triangle( tprlm.p1_, tprlm.p3_, tprlm.p4_, 0, tprlm.multiplicity_);
							continue;
						}

						// not possible: ( (edge1 == e_it_1) and (edge2 == e_it_6) )

						// p3 shared point, p4->p3 shared edge
						if ( (edge1 == e_it_2) and (edge2 == e_it_3) )
						{
							insert_triangle( tprlm.p3_, tprlm.p4_, tprlm.p1_, 0, tprlm.multiplicity_);
							insert_triangle( tprlm.p3_, tprlm.p4_, tprlm.p2_, 0, tprlm.multiplicity_);
							continue;
						}

						// not possible
						// if ( (edge1 == e_it_2) and (edge2 == e_it_4) )

						// p1 shared point, p1->p2 shared edge
						if ( (edge1 == e_it_2) and (edge2 == e_it_5) )
						{
							insert_triangle( tprlm.p1_, tprlm.p2_, tprlm.p3_, 0, tprlm.multiplicity_);
							insert_triangle( tprlm.p1_, tprlm.p2_, tprlm.p4_, 0, tprlm.multiplicity_);
							continue;
						}

						// p3 shared point, p2->p3 shared edge
						if ( (edge1 == e_it_2) and (edge2 == e_it_6) )
						{
							insert_triangle( tprlm.p3_, tprlm.p2_, tprlm.p3_, 0, tprlm.multiplicity_);
							insert_triangle( tprlm.p3_, tprlm.p2_, tprlm.p4_, 0, tprlm.multiplicity_);
							continue;
						}

						// p2 shared point, p2->p1 shared edge
						if ( (edge1 == e_it_3) and (edge2 == e_it_4) )
						{
							insert_triangle( tprlm.p2_, tprlm.p1_, tprlm.p3_, 0, tprlm.multiplicity_);
							insert_triangle( tprlm.p2_, tprlm.p1_, tprlm.p4_, 0, tprlm.multiplicity_);
							continue;
						}

						// not possible: ( (edge1 == e_it_3) and (edge2 == e_it_5) )

						// p3 shared point, p3->p1 shared edge
						if ( (edge1 == e_it_3) and (edge2 == e_it_6) )
						{
							insert_triangle( tprlm.p3_, tprlm.p1_, tprlm.p2_, 0, tprlm.multiplicity_);
							insert_triangle( tprlm.p3_, tprlm.p1_, tprlm.p4_, 0, tprlm.multiplicity_);
							continue;
						}

						// p4 shared point, p4->p3 shared edge
						if ( (edge1 == e_it_4) and (edge2 == e_it_5) )
						{
							insert_triangle( tprlm.p4_, tprlm.p3_, tprlm.p1_, 0, tprlm.multiplicity_);
							insert_triangle( tprlm.p4_, tprlm.p3_, tprlm.p2_, 0, tprlm.multiplicity_);
							continue;
						}

						// p4 shared point, p4->p1 shared edge
						if ( (edge1 == e_it_4) and (edge2 == e_it_6) )
						{
							insert_triangle( tprlm.p4_, tprlm.p1_, tprlm.p3_, 0, tprlm.multiplicity_);
							insert_triangle( tprlm.p4_, tprlm.p1_, tprlm.p2_, 0, tprlm.multiplicity_);
							continue;
						}

						// p4 shared point, p4->p2 shared edge
						if ( (edge1 == e_it_5) and (edge2 == e_it_6) )
						{
							insert_triangle( tprlm.p4_, tprlm.p2_, tprlm.p3_, 0, tprlm.multiplicity_);
							insert_triangle( tprlm.p4_, tprlm.p2_, tprlm.p1_, 0, tprlm.multiplicity_);
							continue;
						}

						throw std::logic_error("Case not caught in rim edge iso-tetrahedron");
					}
				}

				// for less neighbors it gets increasingly difficult to make a sensible choice. we simply add
				// two triangles.
				if ( numNeighbors < 2 )
				{
					insert_triangle( tprlm.p1_, tprlm.p3_, tprlm.p4_, 0, tprlm.multiplicity_);
					insert_triangle( tprlm.p2_, tprlm.p3_, tprlm.p4_, 0, tprlm.multiplicity_);
					continue;
				}

			}

			// here we define a map that turns the point index into the position in kIso_
			// this allows to use strict weak ordering in kIso_ and later redirect the triangles
			int nkIso = bandIsoSurface.size();
			std::vector<int> pointIndexToPosInArray(nkIso);

			kIso_[iw*numBands_+ib] = std::make_shared<std::vector<double>>();
			kIso_[iw*numBands_+ib]->resize(nkIso*3);
			kIsoWeights_[iw*numBands_+ib].resize(nkIso);
			kIsoReducibleWeights_[iw*numBands_+ib].resize(nkIso);
			int c = 0;
			for ( auto const & ngp : bandIsoSurface )
			{
				for ( int i = 0 ; i < 3 ; ++i )
					(*kIso_[iw*numBands_+ib])[c*3+i] = ngp.first.k_[i];
				kIsoWeights_[iw*numBands_+ib][c] = ngp.second.weight;
				kIsoReducibleWeights_[iw*numBands_+ib][c] = ngp.second.reducible_weight;
				pointIndexToPosInArray[ngp.second.pointIndex] = c;
				c++;
			}

			// only now is the order of points clear and
			// in a final step, we define the triangles by linking to
			// the point indices.
			triangles_[iw*numBands_+ib].reserve(triangleToPointIndexMap.size()/3);
			for (int it = 0 ; it < triangleToPointIndexMap.size()/3; ++it)
			{
				int ip1 = pointIndexToPosInArray[triangleToPointIndexMap[it*3+0]];
				int ip2 = pointIndexToPosInArray[triangleToPointIndexMap[it*3+1]];
				int ip3 = pointIndexToPosInArray[triangleToPointIndexMap[it*3+2]];
				triangles_[iw*numBands_+ib].push_back( Triangle(ip1, ip2, ip3, kIso_[iw*numBands_+ib]) );
			}
		}

}

void
TetrahedraIsosurface::compute_iso_vectors(
		LatticeStructure::Tetrahedron const & t,
		std::multimap<double,int> const & ecorner,
		double e,
		std::vector<double> & kIso) const
{
	std::vector<double> p0123;
	t.compute_corner_points(p0123);
	// note: the mapped-to index of ecorner is the index of the corner point

	assert(ecorner.size() == 4);
	auto it = ecorner.begin();
	double e1 = it->first;
	int ip1 = it->second;
	++it;
	double e2 = it->first;
	int ip2 = it->second;
	++it;
	double e3 = it->first;
	int ip3 = it->second;
	++it;
	double e4 = it->first;
	int ip4 = it->second;

	const double eps = 1e-6;

	// triangle
	if ( (e > e1) and (e <= e2) )
	{
		// does not cross e2 <-> e3
		// does not cross e2 <-> e4
		// does not cross e3 <-> e4
		kIso.resize(3*3);
		for ( int i = 0 ; i < 3 ; ++i ) // e1 <-> e2
			kIso[0*3 + i] = std::abs(e2-e1) > eps ?
					p0123[ip1*3+i] + (e-e1)/(e2-e1)*(p0123[ip2*3+i]-p0123[ip1*3+i])
					: 0.5*(p0123[ip1*3+i] + p0123[ip2*3+i]);
		for ( int i = 0 ; i < 3 ; ++i ) // e1 <-> e3
			kIso[1*3 + i] =  std::abs(e3-e1) > eps ?
					p0123[ip1*3+i] + (e-e1)/(e3-e1)*(p0123[ip3*3+i]-p0123[ip1*3+i])
					: 0.5*(p0123[ip1*3+i] + p0123[ip3*3+i]);
		for ( int i = 0 ; i < 3 ; ++i ) // e1 <-> e4
			kIso[2*3 + i] = std::abs(e4-e1) > eps ?
					p0123[ip1*3+i] + (e-e1)/(e4-e1)*(p0123[ip4*3+i]-p0123[ip1*3+i])
					: 0.5*(p0123[ip1*3+i] + p0123[ip4*3+i]);
	}

	// Quadrilateral
	if ( (e > e2) and (e <= e3) )
	{
		// does not cross e3 <-> e4
		// does not cross e1 <-> e2
		kIso.resize(4*3);
		for ( int i = 0 ; i < 3 ; ++i ) // e2 <-> e3
			kIso[0*3 + i] = std::abs(e3-e2) > eps ?
					p0123[ip2*3+i] + (e-e2)/(e3-e2)*(p0123[ip3*3+i]-p0123[ip2*3+i])
					: 0.5*(p0123[ip2*3+i] + p0123[ip3*3+i]);
		for ( int i = 0 ; i < 3 ; ++i ) // e2 <-> e4
			kIso[1*3 + i] = std::abs(e4-e2) > eps ?
					p0123[ip2*3+i] + (e-e2)/(e4-e2)*(p0123[ip4*3+i]-p0123[ip2*3+i])
					: 0.5*(p0123[ip2*3+i] + p0123[ip4*3+i]);
		for ( int i = 0 ; i < 3 ; ++i ) // e1 <-> e4
			kIso[2*3 + i] = std::abs(e4-e1) > eps ?
					p0123[ip1*3+i] + (e-e1)/(e4-e1)*(p0123[ip4*3+i]-p0123[ip1*3+i])
					: 0.5*(p0123[ip4*3+i] + p0123[ip1*3+i]);
		for ( int i = 0 ; i < 3 ; ++i ) // e1 <-> e3
			kIso[3*3 + i] = std::abs(e3-e1) > eps ?
					p0123[ip1*3+i] + (e-e1)/(e3-e1)*(p0123[ip3*3+i]-p0123[ip1*3+i])
					: 0.5*(p0123[ip3*3+i] + p0123[ip1*3+i]);
	}

	// triangle
	if ( (e > e3) and (e <= e4) )
	{
		// does not cross e1 <-> e2
		// does not cross e2 <-> e3
		// does not cross e1 <-> e3
		kIso.resize(3*3);
		for ( int i = 0 ; i < 3 ; ++i ) // e1 <-> e4
			kIso[0*3 + i] = std::abs(e1-e4) > eps ?
					p0123[ip4*3+i] + (e-e4)/(e1-e4)*(p0123[ip1*3+i]-p0123[ip4*3+i])
					: 0.5*(p0123[ip1*3+i] + p0123[ip4*3+i]);
		for ( int i = 0 ; i < 3 ; ++i ) // e2 <-> e4
			kIso[1*3 + i] = std::abs(e2-e4) > eps ?
					p0123[ip4*3+i] + (e-e4)/(e2-e4)*(p0123[ip2*3+i]-p0123[ip4*3+i])
					: 0.5*(p0123[ip2*3+i] + p0123[ip4*3+i]);
		for ( int i = 0 ; i < 3 ; ++i ) // e3 <-> e4
			kIso[2*3 + i] = std::abs(e3-e4) > eps ?
					p0123[ip4*3+i] + (e-e4)/(e3-e4)*(p0123[ip3*3+i]-p0123[ip4*3+i])
					: 0.5*(p0123[ip3*3+i] + p0123[ip4*3+i]);
	}
}

int
TetrahedraIsosurface::get_num_iso_pts(int isoNum, int ibnd) const
{
	assert( (isoNum >= 0) and (isoNum < numIsoE_) );
	assert( (ibnd >= 0) and (ibnd < numBands_) );
	return kIso_[isoNum*numBands_+ibnd]->size()/3;
}


int
TetrahedraIsosurface::get_num_triangles(int isoNum, int ibnd) const
{
	return this->get_triangles(isoNum, ibnd).size();
}

std::vector<TetrahedraIsosurface::Triangle> const &
TetrahedraIsosurface::get_triangles(int isoNum, int ibnd) const
{
	assert( (isoNum >= 0) and (isoNum < numIsoE_) );
	assert( (ibnd >= 0) and (ibnd < numBands_) );
	return triangles_[isoNum*numBands_+ibnd];
}

int
TetrahedraIsosurface::get_nBnd() const
{
	return numBands_;
}

int
TetrahedraIsosurface::get_num_iso_energies() const
{
	return numIsoE_;
}

std::shared_ptr<const LatticeStructure::TetrahedraGrid>
TetrahedraIsosurface::get_tetra_grid() const
{
	return tetra_;
}

void
TetrahedraIsosurface::generate_vtk_polydata(
		int isoNum,
		vtkSmartPointer<vtkPolyData> & polydata) const
{
	auto const & lattice = tetra_->get_grid()->get_lattice();
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
	int idBandOffset = 0;
	for ( int ib = 0 ; ib < numBands_ ; ++ib )
	{
		// Setup points
		for ( int iso = 0 ; iso < this->get_num_iso_pts(isoNum, ib) ; ++iso )
		{
			std::vector<double> coords{ (*kIso_[isoNum*numBands_+ib])[iso*3 + 0],
										(*kIso_[isoNum*numBands_+ib])[iso*3 + 1],
										(*kIso_[isoNum*numBands_+ib])[iso*3 + 2]};
			lattice.reci_direct_to_cartesian(coords);
			points->InsertNextPoint(coords[0], coords[1], coords[2]);
		}

		// Create and insert triangles
		for ( int it = 0 ; it < this->get_num_triangles(isoNum, ib) ; ++it )
		{
			vtkSmartPointer<vtkTriangle> triangle =
					vtkSmartPointer<vtkTriangle>::New();
			auto const & t = triangles_[isoNum*numBands_+ib][it];
			triangle->GetPointIds()->SetId(0, idBandOffset+t.ip1_);
			triangle->GetPointIds()->SetId(1, idBandOffset+t.ip2_);
			triangle->GetPointIds()->SetId(2, idBandOffset+t.ip3_);
			triangles->InsertNextCell(triangle);
		}

		idBandOffset += this->get_num_iso_pts(isoNum, ib);
	}

	// Create a polydata object and add everything to it
	polydata = vtkSmartPointer<vtkPolyData>::New();
	polydata->SetPoints(points);
	polydata->SetPolys(triangles);
}

void
TetrahedraIsosurface::get_irreducible_iso_vector_integration_weights(
		int isoNum,
		int ibnd,
		std::vector<double> & kIso,
		std::vector<double> & kIsoWeights) const
{
	assert( (isoNum >= 0) and (isoNum < numIsoE_) );
	assert( (ibnd >= 0) and (ibnd < numBands_) );
	kIso.assign(kIso_[isoNum*numBands_+ibnd]->begin(), kIso_[isoNum*numBands_+ibnd]->end());
	kIsoWeights.assign(kIsoWeights_[isoNum*numBands_+ibnd].begin(), kIsoWeights_[isoNum*numBands_+ibnd].end());
}

void
TetrahedraIsosurface::get_irreducible_iso_vector_integration_weights_no_multiplicty(
		int isoNum,
		int ibnd,
		std::vector<double> & kIso,
		std::vector<double> & kIsoWeights) const
{
	assert( (isoNum >= 0) and (isoNum < numIsoE_) );
	assert( (ibnd >= 0) and (ibnd < numBands_) );
	kIso.assign(kIso_[isoNum*numBands_+ibnd]->begin(), kIso_[isoNum*numBands_+ibnd]->end());
	kIsoWeights.assign(kIsoReducibleWeights_[isoNum*numBands_+ibnd].begin(),
						kIsoReducibleWeights_[isoNum*numBands_+ibnd].end());
}

void
TetrahedraIsosurface::get_reducible_iso_vector_integration_weights(
		int isoNum,
		int ibnd,
		std::vector<double> & kIso,
		std::vector<double> & kIsoWeights) const
{
	assert(tetra_);
	assert( (isoNum >= 0) and (isoNum < numIsoE_) );
	assert( (ibnd >= 0) and (ibnd < numBands_) );

	int nkIsoIrred = kIso_[isoNum*numBands_+ibnd]->size() / 3;
	assert(kIsoReducibleWeights_[isoNum*numBands_+ibnd].size( ) == nkIsoIrred);
	double eps = tetra_->get_grid()->get_grid_prec();
	auto const & symmertry = tetra_->get_grid()->get_symmetry();

	std::vector<double> allIrredKPts(
			kIso_[isoNum*numBands_+ibnd]->begin(),
			kIso_[isoNum*numBands_+ibnd]->end());

	std::map<NonGridPoint, detail::kinfo> reducibleKPSet;
	for ( int isym = 0 ; isym < symmertry.get_num_symmetries(); ++isym)
	{
		auto reducibleKpts = allIrredKPts;
		symmertry.apply(isym, reducibleKpts.begin(), reducibleKpts.end(), /*lattice periodic = */true);

		for ( int ikIso = 0 ; ikIso < nkIsoIrred ; ++ikIso)
		{
			NonGridPoint ngp( 	reducibleKpts[ikIso*3+0],
								reducibleKpts[ikIso*3+1],
								reducibleKpts[ikIso*3+2],
								eps);
			auto ret = reducibleKPSet.insert( std::make_pair(ngp, detail::kinfo()) );
			ret.first->second.reducible_weight += kIsoReducibleWeights_[isoNum*numBands_+ibnd][ikIso];
		}
	}

	kIso.resize(reducibleKPSet.size()*3);
	kIsoWeights.resize(reducibleKPSet.size());
	int ikRed = 0;
	for (auto const & p : reducibleKPSet )
	{
		for (int i = 0 ; i < 3 ; ++i)
			kIso[ikRed*3+i] = p.first.k_[i];
		kIsoWeights[ikRed] = p.second.reducible_weight;
		++ikRed;
	}
}

void
TetrahedraIsosurface::write_polydata_file( std::string const & filename,
		int isoNum,
		std::shared_ptr<std::vector<double>> data) const
{
	vtkSmartPointer<vtkPolyData> polydata;
	this->generate_vtk_polydata(isoNum, polydata);

	vtkSmartPointer<vtkXMLPolyDataWriter> xmlWriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	xmlWriter->SetInputData(polydata);
	xmlWriter->SetFileName( filename.c_str() );
	xmlWriter->Update();
}

namespace detail
{
bool operator< (NonGridPoint const & ngp1, NonGridPoint const & ngp2)
{
	assert(ngp1.eps_ == ngp2.eps_);
	for ( int i = 0 ; i < 3 ; ++i)
		if ( std::fabs(ngp1.k_[i] - ngp2.k_[i]) > ngp1.eps_ )
			return ngp1.k_[i] < ngp2.k_[i];
	return false;
};


double
Triangle::area() const
{
	assert(pts_);
	return Algorithms::helperfunctions::triangle_area(
			this->get_point_coords(0),
			this->get_point_coords(1),
			this->get_point_coords(2));
}

std::vector<double>
Triangle::get_point_coords(int i) const
{
	assert(pts_);
	assert( (i >= 0) and (i < 3));
	int id = ip1_;
	if ( i == 1 )
		id = ip2_;
	if ( i == 2 )
		id = ip3_;
	assert( (id*3 >= 0) and (pts_->size() >= id*3+3));
	return std::vector<double>( &((*pts_)[id*3+0]), &((*pts_)[id*3+0]) + 3);
}

}

} /* namespace ElectronicStructure */
} /* namespace elephon */
