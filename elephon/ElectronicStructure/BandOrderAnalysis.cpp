/*	This file BandOrderAnalysis.cpp is part of elephon.
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
 *  Created on: Sep 11, 2017
 *      Author: A. Linscheid
 */

#include "ElectronicStructure/BandOrderAnalysis.h"
#include "Algorithms/FFTInterface.h"
#include "IOMethods/WriteVASPRealSpaceData.h"
#include <map>
#include <set>
#include <iostream>

namespace elephon
{
namespace ElectronicStructure
{

void
BandOrderAnalysis::compute_band_order_overlap(
		ElectronicBands const & bands,
		Wavefunctions const & wfct)
{
	nB_ = wfct.get_num_bands();
	auto & kgrid = wfct.get_k_grid();
	auto d = kgrid.get_grid_dim();
	bandOrder_.assign( nB_*kgrid.get_np_irred(), -1 );

	// A note on the Algorithm: the idea is to propagate the ordering
	// of a seed (e.g. the first k point) through the mesh.
	// Problem:
	// Degeneracies have to be dealt with, since for such bands the eigenvector
	// can be any linear combination, thus a overlap ordering is ill defined.
	// Solution:
	// The current implementation checks for degenerate eigenvalues at the seed k point
	// (the center of a 3x3x3 subset in the regular mesh) and for every degenerate band
	// attempts to find a neighbor with confirmed order where it is not degenerate that replaces
	// the seed. This should be able to deal with degenerate points and lines. Degenerate planes
	// or even volumes would require a non-local probe.

	// lambda to check if bands are degenerate.
	const double degeneracy_thrhld = 1e-3;
	auto check_degenerate_bands = [&] (int ikir, int ib){
		// Since bands by default are ordered in energy, we only check neighbors
		if ( ib + 1 < nB_ )
			if ( std::abs(bands(ikir, ib) - bands(ikir, ib+1)) < degeneracy_thrhld )
				return true;
		if ( ib - 1 >= 0 )
			if ( std::abs(bands(ikir, ib) - bands(ikir, ib-1)) < degeneracy_thrhld )
				return true;
		return false;
	};

	auto fftMax = wfct.get_max_fft_dims();

	// step 1: set up a structure which tells for a k point all its neighbors up to 2 away in the grid
	// 			each element of the set points to a index of the map and thereby establishes neighbors.
	// choose the irreducible wedge with the smallest x, then smallest y then smallest z.
	// because the k grid is ordered x major, this means we simply search the smallest ikr in the star of k
	auto irreducible_wedge = [&kgrid] (int ikir){
		std::pair<int, int> symLow = std::make_pair(0, std::numeric_limits<int>::max() );
		for ( int istar = 0 ; istar < kgrid.get_maps_sym_irred_to_reducible()[ikir].size(); ++istar)
		{
			int ikrTest = kgrid.get_maps_irreducible_to_reducible()[ikir][istar];
			if ( ikrTest < symLow.second )
				symLow = std::make_pair(istar, ikrTest);
		}
		return kgrid.get_maps_irreducible_to_reducible()[ikir][symLow.first];
	};
	std::map<int, std::vector<int>> linkedGrid;
	for ( int ikir = 0; ikir < kgrid.get_np_irred(); ++ikir )
	{
		int ikr = irreducible_wedge(ikir);
		auto xyz = kgrid.get_reducible_to_xyz(ikr);
		auto xyzMod = xyz;
		std::map<int, int> neighbors;
		for ( int pz = -2; pz <= 2; ++pz)
			for ( int py = -2; py <= 2; ++py)
				for ( int px = -2; px <= 2; ++px)
				{
					if ( (px == 0) and (py == 0) and (pz == 0) )
						continue;
					xyzMod[0] = xyz[0]+px;
					xyzMod[1] = xyz[1]+py;
					xyzMod[2] = xyz[2]+pz;
					for ( int i = 0 ; i < 3 ; ++i )
					{
						xyzMod[i] = xyzMod[i] < d[i] ? xyzMod[i] : xyzMod[i] - d[i];
						xyzMod[i] = xyzMod[i] >= 0 ? xyzMod[i] : xyzMod[i] + d[i];
					}
					int neighIkR = kgrid.get_xyz_to_reducible(xyzMod);
					int neigIkIr = kgrid.get_maps_red_to_irreducible()[neighIkR];
					// remove the point itself from its neighbors even for periodic boundary conditions, symmetry etc
					if ( ikr == irreducible_wedge(neigIkIr) )
						continue;
					int distanceOrder = abs(px)+abs(py)+abs(pz);
					auto el = neighbors.find(neigIkIr);
					neighbors[neigIkIr] = (el == neighbors.end()) ? distanceOrder : std::min(neighbors[neigIkIr], distanceOrder);
				}

		std::multimap<int, int> uniqNeigDistanceSorted;
		for ( auto n : neighbors)
			uniqNeigDistanceSorted.insert(std::make_pair(n.second, n.first));

		// we map from irreducible back to reducible. Thus uniqNeighborsDistanceSorted contains reducible
		// indices of the irreducible wedge of the BZ. The must be by logic valid keys of the map once we have completed this
		// step.
		std::vector<int> uniqNeighborsDistanceSorted;
		uniqNeighborsDistanceSorted.reserve(uniqNeighborsDistanceSorted.size());
		for ( auto n : uniqNeigDistanceSorted )
			uniqNeighborsDistanceSorted.push_back(irreducible_wedge(n.second));
		linkedGrid[ikr] = std::move(uniqNeighborsDistanceSorted);
	}

	// check the logic that neighbors are also in the irredcuble wedge
	for ( auto k : linkedGrid )
		for ( auto n : k.second )
			if ( linkedGrid.find(n) == linkedGrid.end() )
				throw std::logic_error("Linked grid neighbors lead out of the irreducible zone");
	if ( linkedGrid.size() != kgrid.get_np_irred() )
		throw std::logic_error("Linked grid does not span the irreducible zone");

	// Build a set of paths of points. The first one starts at ikr = 0.
	std::vector<std::pair<int, std::vector<int>>> paths;
	bool lastIn = false;
	for ( int ikr = 0; ikr < kgrid.get_np_red(); ++ikr )
	{
		if ( linkedGrid.find(ikr) == linkedGrid.end() )
		{
			lastIn = false;
			continue; // not in the irreducible zone
		}
		if (not lastIn)
		{
			paths.push_back(std::make_pair(ikr, std::vector<int>{ikr}));
			lastIn = true;
			continue;
		}
		paths.rbegin()->second.push_back(ikr);
	}

	// place a seed order
	std::vector<bool> bandConventionFixed(nB_, false);
	for ( int ib = 0 ; ib < nB_; ++ib)
	{
		bandOrder_[ib] = ib;
		if ( not check_degenerate_bands(0, ib) )
			bandConventionFixed[ib] = true;
	}

	const double significant_overlap = 0.05;

	// step 2: walk through the structure and compare with neighbors.
	//	if a point has degenerate bands, replace the wavefunction with one of the closest neighbors for these bands
	std::vector<std::complex<float>> wfct1c(fftMax[0]*fftMax[1]*fftMax[2]);
	double progress = 0.0;
	for ( auto kpath : paths )
	{
		std::cout << "Progress: "<< progress*100 << "%"<< std::endl;
		progress += double(kpath.second.size())/kgrid.get_np_irred();
		for ( int ikr : kpath.second )
		{
			auto it = linkedGrid.find(ikr);
			auto const & neighbors = it->second;
			int ikir = kgrid.get_maps_red_to_irreducible()[ikr];

			// in the set of this point and its neighbors classify into
			// points already set and those that are not. We only look at band 0.
			std::vector<int> sources;
			std::vector<int> targets;
			if ( bandOrder_[ikir*nB_] >= 0 )
				sources.push_back(ikr);
			else
				targets.push_back(ikr);
			for (int n : neighbors )
				if ( bandOrder_[kgrid.get_maps_red_to_irreducible()[n]*nB_] >= 0 )
					sources.push_back(n);
				else
					targets.push_back(n);

			assert(sources.size() > 0);
			if ( targets.size() == 0 )
				continue;

			std::vector<std::map<float,int>> overlapMatrixNeigbors(nB_*targets.size());
			std::vector<std::pair<int,float>> bandsEnergySource(nB_);
			for ( int ib = 0 ; ib < nB_; ++ib)
			{
				auto s = sources.begin();
				while( check_degenerate_bands(kgrid.get_maps_red_to_irreducible()[*s], ib) )
				{
					if ( (++s) == sources.end() )
					{
						if ( bandConventionFixed[ib] )
							throw std::logic_error(std::string("Unable to find a source for "
									"propagating a band convention for band ")+std::to_string(ib));
						// here we treat the (may be not so) special case that the first k points have degenerate bands
						// in that case, we promote a target with non-degenerate bands to the source and mark the convention as fixed.
						auto t = targets.begin();
						while( check_degenerate_bands(kgrid.get_maps_red_to_irreducible()[*t], ib) )
						{
							if ( (++t) == targets.end() )
								throw std::runtime_error(std::string("Unable to fix band convention"
										" for band ")+std::to_string(ib)+". This can happen when a band is highly degenerate."
												" Try different grids.");
						}
						bandOrder_[(kgrid.get_maps_red_to_irreducible()[*t])*nB_+ib]
								   = bandOrder_[(kgrid.get_maps_red_to_irreducible()[*s])*nB_+ib];
						bandConventionFixed[ib] = true;
						s = t;
						break;
					}
				}
				int iksr = *s;

				// first load the source wfct and save it in a regular mesh
				std::vector<std::vector<int>> fftmaps;
				auto ks = kgrid.get_vector_direct(iksr);
				wfct.compute_Fourier_maps(ks, fftmaps);
				int npws = fftmaps[0].size()/3;
				auto fms = fftmaps[0];

				std::vector<std::vector<std::complex<float>>> reducWfcts;
				std::vector<int> npwPerK;
				wfct.generate_reducible_grid_wfcts(std::vector<int>{ib}, std::vector<int>{iksr}, reducWfcts, npwPerK);
				assert( npws == npwPerK[0] );
				std::fill(wfct1c.begin(), wfct1c.end(), std::complex<float>(0.0f));
				for (int ig = 0; ig < npws; ++ig)
					wfct1c[fms[ig*3+0]+fftMax[0]*(fms[ig*3+1]+fftMax[1]*fms[ig*3+2])] =
							std::conj(reducWfcts[0][ig]);

				int iksir = kgrid.get_maps_red_to_irreducible()[iksr];
				bandsEnergySource[ib] = std::make_pair(iksir, bands(iksir, ib));

				for ( auto itt = targets.begin(); itt != targets.end(); ++itt)
				{
					int iktr = *itt;
					if ( iksr == iktr )
						continue;

					auto kt = kgrid.get_vector_direct(iktr);
					wfct.compute_Fourier_maps(kt, fftmaps);
					int npwt = fftmaps[0].size()/3;
					auto const & fmt = fftmaps[0];

					int nei_d = std::distance(targets.begin(), itt);


					// get a list with bands sorted by closest to farthest in the current order
					// the hope is that matching bands are close in energy, too.
					std::vector<int> bandIndices(nB_);
					for ( int ibp = 0 ; ibp < nB_; ++ibp)
						bandIndices[ibp] = ibp;
					std::vector<int> lbdist = bandIndices;
					for ( auto &b : lbdist)
						b -= ib;
					std::sort(lbdist.begin(), lbdist.end(), [](int a, int b){return abs(a)<abs(b); } );

					// here comes the part where we propagate the band convention from the k point iksr
					// to the target k point iktr
					std::map<float,int> ovlap;
					float sum = 0.0f;
					for ( int ibo : lbdist)
					{
						int ibp = ib + ibo;
						assert((ibp >= 0) and (ibp < nB_));
						wfct.generate_reducible_grid_wfcts(std::vector<int>{ibp}, std::vector<int>{iktr}, reducWfcts, npwPerK);
						assert( npwt == npwPerK[0] );
						std::complex<float> overlap(0.0f);
						for (int ig = 0; ig < npwt; ++ig)
							overlap += wfct1c[fmt[ig*3+0]+fftMax[0]*(fmt[ig*3+1]+fftMax[1]*fmt[ig*3+2])]
											 *reducWfcts[0][ig];
						sum += std::abs(overlap);
						ovlap.insert(std::make_pair(std::abs(overlap), ibp));
						// Unfortunately sum is not normalized to 1 because the wavefunctions are not
						// normalized. However, they are close ~5% so this should work in practice.
//						if ( ovlap.rbegin()->first > (1.0f-sum) )
//							break;
					}
					int iktir = kgrid.get_maps_red_to_irreducible()[iktr];
					// handle the case where the band in not degenerate, i.e. no problem
					if ( (ovlap.rbegin()->first > significant_overlap)
							and not check_degenerate_bands(iktir, ovlap.rbegin()->second))
						bandOrder_[iktir*nB_+ovlap.rbegin()->second] = bandOrder_[iksir*nB_+ib];
					// for high energy bands, it is possible that they are not continued and another one
					// takes over. Also if the band overlaps into a degenerate set, we simply order according to heuristics
					overlapMatrixNeigbors[ib+nB_*nei_d] = std::move(ovlap);
				}
			}

			// at this point we need to fill gaps where target bands are degenerate or bands leave the
			// energy window and another band takes over.
			for ( auto itt = targets.begin(); itt != targets.end(); ++itt)
			{
				int iktr = *itt;
				int iktir = kgrid.get_maps_red_to_irreducible()[iktr];

				int nei_d = std::distance(targets.begin(), itt);

				// find degenerate subspaces and store the band indices in each one
				std::vector<std::vector<int>> degenSubsets;
				for ( int ib = 0 ; ib < nB_; ++ib)
					if ( ib + 1 < nB_ )
						if ( std::abs(bands(ikir, ib) - bands(ikir, ib+1)) < degeneracy_thrhld )
						{
							if ( degenSubsets.empty() )
								degenSubsets.push_back( std::vector<int>{ib, ib+1} );
							else
							{
								if ( *degenSubsets.rbegin()->rbegin() == ib )
									degenSubsets.rbegin()->push_back(ib+1);
								else
									degenSubsets.push_back( std::vector<int>{ib, ib+1} );
							}
						}

				// consult the overlap matrix if bands from the source overlap
				// with a given subspace. They are assigned in ascending order
				for ( auto & deg : degenSubsets)
				{
					std::set<int> sourceBandsThatOverlap;
					for ( int ib : deg)
						for ( int ibp = 0 ; ibp < nB_; ++ibp)
						{
							auto itr = overlapMatrixNeigbors[ibp+nB_*nei_d].begin();
							for ( ; itr !=  overlapMatrixNeigbors[ibp+nB_*nei_d].end(); ++itr)
							{
								if ( (itr->second == ib) and (itr->first >= significant_overlap) )
								{
									// significant overlap into degenerate subspace
									sourceBandsThatOverlap.insert(ibp);
								}
							}
						}
					auto it = sourceBandsThatOverlap.begin();
					for ( int ib : deg)
					{
						if ( it == sourceBandsThatOverlap.end())
							break;
						bandOrder_[iktir*nB_+ib] = *(it++);
					}
				}

				// last resort: use the closest in energy band of the source that has not been assigned.
				std::set<int> assigned;
				for ( int ib = 0 ; ib < nB_; ++ib)
					if ( bandOrder_[iktir*nB_+ib] >= 0 )
						assigned.insert(bandOrder_[iktir*nB_+ib]);

				for ( int ib = 0 ; ib < nB_ ; ++ib )
				{
					if ( bandOrder_[iktir*nB_+ib] < 0 )
					{
						std::map<float, int> closestUnassignedBands;
						for ( int ibp = 0 ; ibp < nB_ ; ++ibp )
							if ( assigned.find(ibp) == assigned.end() )
								closestUnassignedBands.insert(
									std::make_pair(std::abs(bands(iktir, ib)-bandsEnergySource[ibp].second), ibp));
						assert(closestUnassignedBands.size() > 0);
						bandOrder_[iktir*nB_+ib] = closestUnassignedBands.begin()->second;
					}
				}

				std::vector<int> unassigned;
				unassigned.reserve(nB_ - int(assigned.size()));
				for ( int ib = 0 ; ib < nB_; ++ib)
					if ( assigned.find(ib) == assigned.end() )
						unassigned.push_back(ib);
			}

		}
	}
}

void
BandOrderAnalysis::compute_band_overlap_matrix(
		Wavefunctions const & wfct,
		int ikr1, std::vector<int> const & bands1,
		int ikr2, std::vector<int> const & bands2,
		std::vector<float> & overlaps) const
{
	auto const & kgrid = wfct.get_k_grid();

	std::vector<std::vector<int>> fftmaps;
	auto ks = kgrid.get_vector_direct(ikr1);
	wfct.compute_Fourier_maps(ks, fftmaps);
	int npws = fftmaps[0].size()/3;
	auto fms = fftmaps[0];

	Algorithms::FFTInterface fft;

	auto fftMax = wfct.get_max_fft_dims();
	std::vector<std::complex<float>> wfct1c, wfct2;

	std::vector<int> star1;
	int ikir1 = kgrid.get_maps_red_to_irreducible()[ikr1];
	for (int is1 = 0 ;is1 < kgrid.get_maps_sym_irred_to_reducible()[ikir1].size() ; ++is1)
		star1.push_back(kgrid.get_maps_irreducible_to_reducible()[ikir1][is1]);

	std::vector<int> star2;
	int ikir2 = kgrid.get_maps_red_to_irreducible()[ikr2];
	for (int is2 = 0 ;is2 < kgrid.get_maps_sym_irred_to_reducible()[ikir2].size() ; ++is2)
		star2.push_back(kgrid.get_maps_irreducible_to_reducible()[ikir2][is2]);

	std::vector<std::vector<std::complex<float>>> reducWfcts;
	std::vector<int> npwPerK;
	wfct.generate_reducible_grid_wfcts(bands1, star1, reducWfcts, npwPerK);
	std::fill(wfct1c.begin(), wfct1c.end(), std::complex<float>(0.0f));
	for ( int is =1 ; is < star1.size(); ++is)
	{
		assert(reducWfcts[0].size() == reducWfcts[is].size());
		for ( int i = 0 ; i < reducWfcts[0].size(); ++i )
			reducWfcts[0][i] += reducWfcts[is][i];
	}
	for ( auto & w: reducWfcts[0] )
		w = std::conj(w)/float(star1.size());
	fft.fft_sparse_data(fms, fftMax, reducWfcts[0], nB_, +1, wfct1c, fftMax);

	auto kt = kgrid.get_vector_direct(ikr2);
	wfct.compute_Fourier_maps(kt, fftmaps);
	int npwt = fftmaps[0].size()/3;
	auto const & fmt = fftmaps[0];


	wfct.generate_reducible_grid_wfcts(bands2, star2, reducWfcts, npwPerK);
	for ( int is =1 ; is < star2.size(); ++is)
	{
		assert(reducWfcts[0].size() == reducWfcts[is].size());
		for ( int i = 0 ; i < reducWfcts[0].size(); ++i )
			reducWfcts[0][i] += reducWfcts[is][i];
	}
	for ( auto & w: reducWfcts[0] )
		w /= star2.size();
	fft.fft_sparse_data(fmt, fftMax, reducWfcts[0], nB_, -1, wfct2, fftMax);

	int npts = fftMax[0]*fftMax[1]*fftMax[2];
	overlaps.assign(nB_*nB_, 0.0f);
	for ( int ib1 = 0 ; ib1 < bands1.size(); ++ib1 )
		for (int ib2 = 0 ; ib2 < bands2.size(); ++ib2 )
		{
			if ( overlaps[bands1[ib1]+nB_*bands2[ib2]] != 0.0f )
				continue;

			std::complex<float> c(0.0f);
			auto it1 = wfct1c.begin()+npts*ib1;
			auto it2 = wfct2.begin()+npts*ib2;
			for (; it1 != (wfct1c.begin()+npts*(ib1+1)); ++it1, ++it2)
			{
				c += (*it1)*(*it2);
			}
//			for (int ig = 0; ig < npwt; ++ig)
//				c += wfct1c[fmt[ig*3+0]+fftMax[0]*(fmt[ig*3+1]+fftMax[1]*(fmt[ig*3+2]+fftMax[2]*ib1))]
//								 *reducWfcts[0][ig+npwt*ib2];
			overlaps[bands1[ib1]+nB_*bands2[ib2]] = std::abs(c);
			overlaps[bands2[ib2]+nB_*bands1[ib1]] = overlaps[bands1[ib1]+nB_*bands2[ib2]];
		}

	LatticeStructure::UnitCell uc;
	LatticeStructure::Atom al("Al",{0.0,0.0,0.0},{false,false,false});
	std::vector<LatticeStructure::Atom> atoms(1, al);
	auto sym = wfct.get_k_grid().get_symmetry();
	sym.set_reciprocal_space_sym(false);
	uc.initialize(atoms, wfct.get_k_grid().get_lattice(), sym);

	std::vector<double> data(npts);
	for (int ir = 0 ; ir < npts; ++ir)
	{
		data[ir] = std::real(wfct1c[ir + npts*(nB_-1)]*std::conj(wfct1c[ir + npts*(nB_-1)]));
	}

	IOMethods::WriteVASPRealSpaceData writer;
	writer.write_file("/tmp/LOCPOT", "bla", fftMax, uc, data);

	for (int ir = 0 ; ir < npts; ++ir)
	{
		data[ir] = std::real(wfct2[ir + npts*(nB_-1)]*std::conj(wfct2[ir + npts*(nB_-1)]));
	}
\
	writer.write_file("/tmp/CHGCAR", "bla", fftMax, uc, data);
}

int
BandOrderAnalysis::operator() (int ikir, int ib) const
{
	assert((ikir*nB_+ib >= 0) and (ikir*nB_+ib<bandOrder_.size()));
	return bandOrder_[ikir*nB_+ib];
}

} /* namespace ElectronicStructure */
} /* namespace elephon */
