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
#include <stdexcept>

namespace elephon
{
namespace ElectronicStructure
{

void
BandOrderAnalysis::compute_band_order_overlap(
		ElectronicBands const & bands,
		Wavefunctions const & wfct)
{
	throw std::logic_error("not implemented. The overlap algorithm appears to be a weak classification.");

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


	auto fftMax = wfct.get_max_fft_dims();

	std::map<int, std::pair<bool, std::vector<int> >> linkedGrid;
	for ( int ikir = 0; ikir < kgrid.get_np_irred(); ++ikir )
	{
		int ikr = kgrid.get_maps_irreducible_to_reducible()[ikir][0];
		auto xyz = kgrid.get_reducible_to_xyz(ikr);
		auto xyzMod = xyz;
		std::map<int, int> neighbors;
		for ( int pz = -1; pz <= 1; ++pz)
			for ( int py = -1; py <= 1; ++py)
				for ( int px = -1; px <= 1; ++px)
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
					if ( ikr == kgrid.get_maps_irreducible_to_reducible()[neigIkIr][0] )
						continue;
					auto xyz2 = kgrid.get_reducible_to_xyz(kgrid.get_maps_irreducible_to_reducible()[neigIkIr][0]);
					int distanceOrder = 0;
					for (int i = 0 ; i < 3 ; ++i)
					{
						distanceOrder += abs(xyz2[i]-xyz[i]);
					}
					if ( distanceOrder > 3 )
						continue;
					auto el = neighbors.find(neigIkIr);
					neighbors[neigIkIr] = (el == neighbors.end()) ? distanceOrder : std::min(neighbors[neigIkIr], distanceOrder);
				}

		std::multimap<int, int> uniqNeigDistanceSorted;
		for ( auto n : neighbors)
			uniqNeigDistanceSorted.insert(std::make_pair(n.second, n.first));

		std::vector<int> uniqNeighborsDistanceSorted;
		uniqNeighborsDistanceSorted.reserve(uniqNeighborsDistanceSorted.size());
		for ( auto n : uniqNeigDistanceSorted )
			uniqNeighborsDistanceSorted.push_back(kgrid.get_maps_irreducible_to_reducible()[n.second][0]);

		bool isHighSymmertryK;
		auto ksym = kgrid.get_symmetry();
		auto k = kgrid.get_vector_direct(ikr);
		ksym.small_group(k);
		isHighSymmertryK = (ksym.get_num_symmetries_no_T_rev() > 1);

		linkedGrid[ikr] = std::move(std::make_pair(isHighSymmertryK, std::move(uniqNeighborsDistanceSorted)) );
	}

	// check the logic that neighbors are also in the irredcuble wedge
	for ( auto k : linkedGrid )
		for ( auto n : k.second.second )
			if ( linkedGrid.find(n) == linkedGrid.end() )
				throw std::logic_error("Linked grid neighbors lead out of the irreducible zone");
	if ( linkedGrid.size() != kgrid.get_np_irred() )
		throw std::logic_error("Linked grid does not span the irreducible zone");

	// Build a set of paths of non-high-symmetry k points.
	std::vector<std::vector<int>> paths;
	bool lastIn = false;
	std::map<int, std::pair<int,int>> assignedKptsPaths;
	bool seed = true;

	auto start_new_path = [&] (int ikr, std::vector<int> const & neighbors)
		{
			int ikn = -1;
			for (  int ik ; ik < neighbors.size() ; ++ik )
			{
				auto test = assignedKptsPaths.find(neighbors[ik]);
				if ( test != assignedKptsPaths.end())
				{
					ikn = neighbors[ik];
					break;
				}
			}
			if ( ikn < 0 )
			{
				if (seed == false )
					throw std::logic_error("Unable to find a previous k point in close proximity");
				paths.push_back(std::vector<int>{ikr});
				seed = false;
			}
			else
			{
				paths.push_back(std::vector<int>{ikn,ikr});
			}

		};

	for ( int ikr = 0; ikr < kgrid.get_np_red(); ++ikr )
	{
		auto el = linkedGrid.find(ikr);
		if ( el == linkedGrid.end() )
		{
			lastIn = false;
			continue; // not in the irreducible zone
		}
		if ( el->second.first )
		{
			lastIn = false;
			continue; // is a high symmetry point
		}

		if (not lastIn) // start a new path
		{
			assignedKptsPaths[ikr] = std::make_pair(int(paths.size()), seed ? 0 : 1);
			start_new_path( ikr, el->second.second );
			lastIn = true;
			continue;
		}
		// check if the last entry is adjacent to the current location in the cell [0,1[
		auto xyzLast = kgrid.get_reducible_to_xyz(*paths.rbegin()->rbegin());
		auto xyzNew = kgrid.get_reducible_to_xyz(ikr);
		double sum = 0;
		for (int i = 0 ; i < 3 ; ++i)
			sum += std::pow(std::abs(xyzNew[i]-xyzLast[i]),2);
		bool close = std::sqrt(sum) <= std::sqrt(3.1);

		if ( not close )
		{
			start_new_path( ikr, el->second.second );
		}

		assignedKptsPaths[ikr] = std::make_pair(int(paths.size())-1, int(paths.rbegin()->size()));
		paths.rbegin()->push_back(ikr);
	}

	std::vector<int> bandorder= {0, 1, 2, 3, 4};
	auto bandorderNext = bandorder;
	std::map<int,std::vector<int>> bandOrderMap;
	bandOrderMap[paths[0][0]] = bandorder;

	std::vector<float> overlaps;
	this->compute_band_overlap_matrix(wfct, 1066, bandorder , 1068, bandorder,overlaps );
	for (int ib  : bandorder )
	{
		for (int ibp  : bandorder )
			std::cout << '\t'<< overlaps[ib*nB_ + ibp] ;
		std::cout <<std::endl;
	}

for (int ip = 0 ; ip < paths.size(); ++ip )
{
	bandorder = bandOrderMap[paths[ip][0]];
	assert(not bandorder.empty());
for (int ikp = 1 ; ikp < paths[ip].size(); ++ikp )
{
	int ikr1 = paths[ip][ikp-1];
	int ikr2 = paths[ip][ikp];
	std::vector<int> b = {0, 1, 2, 3, 4};
	this->compute_band_overlap_matrix(wfct, ikr1, b , ikr2, b,overlaps );
	auto xyz = kgrid.get_reducible_to_xyz(ikr1);
	auto k1 = kgrid.get_vector_direct(ikr1);
	auto k2 = kgrid.get_vector_direct(ikr2);
	double norm = std::sqrt(std::pow(k1[0]-k2[0],2)+std::pow(k1[1]-k2[1],2)+std::pow(k1[2]-k2[2],2));
	int ikir = kgrid.get_maps_red_to_irreducible()[ikr1];
	std::cout << ip << '\t' << ikr1 << '\t' << ikr2<<" **********" << std::endl;
	std::cout << xyz[0] << '\t'<< xyz[1] << '\t'<< xyz[2]<<'\t'<<norm<<'\t';
	std::cout << std::endl;

	std::vector<int> bandMatrix(nB_*nB_, 0);
	bool AssignComplete = true;
	for (int ib  : b )
	{
		bool failAssign = true;
		std::map<float, int> ovlp;
		std::pair<int, float> mp = std::make_pair(ib, 0);
		for (int ibp  : b )
		{
			ovlp.insert(std::make_pair(overlaps[ib*nB_ +ibp],ibp));
			std::cout << '\t'<< overlaps[ib*nB_ +ibp] ;
			if ( (overlaps[ib*nB_ +ibp] > mp.second) and ( overlaps[ib*nB_ +ibp] > 0.8) )
			{
				failAssign = false;
				mp = std::make_pair(ibp, overlaps[ib*nB_ +ibp]);
			}
		}
		if ( not failAssign )
			bandMatrix[ib*nB_ + ovlp.rbegin()->second] = 1;
		AssignComplete = AssignComplete and (not failAssign);
		bandorderNext[ib] = bandorder[mp.first];
		std::cout << std::endl;
	}
	std::vector<int> clmnSm(nB_, 0);
	std::vector<int> rowSm(nB_, 0);
	for (int ib  : b )
	{
		for (int ibp  : b )
		{
			clmnSm[ib] += bandMatrix[ib*nB_ + ibp];
			rowSm[ib] += bandMatrix[ibp*nB_ + ib];
		}
	}
	std::cout << std::boolalpha << AssignComplete <<std::endl;

	for (int ib  : b )
		std::cout << bands(ikir, bandorder[ib]) << '\t';
	std::cout << std::endl;
	for (int ib  : b )
	{
		for (int ibp  : b )
			std::cout << '\t'<< bandMatrix[ib*nB_ + ibp] ;
		std::cout <<std::endl;
	}
	bandorder = bandorderNext;
	bandOrderMap[paths[ip][ikp]] = bandorder;
	std::cout << std::endl;
	std::cout << "**********" << std::endl<< std::endl;
}
}
std::exit(0);
	// step 2: walk through the structure and compare with neighbors.
	//	if a point has degenerate bands, replace the wavefunction with one of the closest neighbors for these bands
	std::vector<std::complex<float>> wfct1c(fftMax[0]*fftMax[1]*fftMax[2]);
	for ( int ib = 0 ; ib < nB_; ++ib)
	{
		bool bandConventionFixed = false;
		std::cout << "Progress: "<< ib*100.0/nB_ << "%"<< std::endl;
		for ( int ip = 0; ip < paths.size(); ++ip )
		{
			for ( int ikr : paths[ip] )
			{
auto xyz = kgrid.get_reducible_to_xyz(ikr);
//std::cout << ip << '\t' << xyz[0]  << '\t' << xyz[1]  << '\t'<<xyz[2] << std::endl;
if ( (xyz[1] != 2) or (xyz[2] != 1 ) )
	continue;
				auto it = linkedGrid.find(ikr);
				auto neighbors = it->second.second;
				int ikir = kgrid.get_maps_red_to_irreducible()[ikr];

				bool high_symmetry_k =
						kgrid.get_maps_sym_irred_to_reducible()[ikir].size() != kgrid.get_symmetry().get_num_symmetries();
				if ( high_symmetry_k )
					continue;

				std::vector<std::pair<int,int>> sources; // kpt, pos of the band index ib
				std::vector<int> targets;
				int ibMap = ib;
				if ( not bandConventionFixed )
				{
					sources.push_back(std::make_pair(ikr,ib));
					for (int n : neighbors )
						targets.push_back(n);
				}
				else
				{
					auto neighborsp = neighbors;
					neighborsp.push_back(ikr);
					for (int n : neighborsp )
					{
auto xyz = kgrid.get_reducible_to_xyz(n);
if ( (xyz[1] != 1) or (xyz[2] != 1 ) )
	continue;
						bool isSource = false;
						for ( int ibp = 0 ; ibp < nB_; ++ibp)
							if ( bandOrder_[kgrid.get_maps_red_to_irreducible()[n]*nB_ + ibp] == ib )
							{
								// make sure high symmetry points are never a source
								int iksir = kgrid.get_maps_red_to_irreducible()[n];
								if ( kgrid.get_maps_sym_irred_to_reducible()[iksir].size() !=
										kgrid.get_symmetry().get_num_symmetries() )
									continue;
								sources.push_back(std::make_pair(n, ibp));
								isSource = true;
								break;
							}
						if (not isSource)
						{
//							auto xyz = kgrid.get_reducible_to_xyz(n);
//							auto xyz0 = kgrid.get_reducible_to_xyz(ikr);
//							if ( (abs(xyz[0]-xyz0[0]) <= 1) and (abs(xyz[1]-xyz0[1]) <= 1) and (abs(xyz[2]-xyz0[2]) <= 1))
								targets.push_back(n);
						}
					}
				}

				if ( sources.size() < 1 )
					throw std::runtime_error("Continuation of bands failed - no source point in range");
				if ( targets.size() == 0 )
					continue;

				int iksr = sources.begin()->first;
				int iksir = kgrid.get_maps_red_to_irreducible()[iksr];

				// first load the source wfct and save it in a regular mesh
				std::vector<std::vector<int>> fftmaps;
				auto ks = kgrid.get_vector_direct(iksr);
				wfct.compute_Fourier_maps(ks, fftmaps);
				int npws = fftmaps[0].size()/3;
				auto fms = fftmaps[0];

				double bnd1 = bands(iksir, ibMap);
				auto k1 = kgrid.get_vector_direct(iksr);
				std::vector<std::vector<int>> candidateCrossings(targets.size());
				for ( int it = 0;  it < targets.size(); ++it)
				{
					int iktr = targets[it];
					if ( iksr == iktr )
						continue;
					int iktir = kgrid.get_maps_red_to_irreducible()[iktr];
					for ( int ibp = 0 ; ibp < nB_; ++ibp)
					{
						// estimate gradient and discard bands unrealistically far away
						double bnd2 = bands(iktir, ibp);
						auto k2 = kgrid.get_vector_direct(iktr);
						for ( int i = 0 ; i < 3 ; ++i)
							k2[i] -= k1[i];
						kgrid.get_lattice().reci_direct_to_cartesian_2pibya(k2);
						double grad = (bnd2-bnd1)/std::sqrt(std::pow(k2[0],2)+std::pow(k2[1],2)+std::pow(k2[2],2));
						if( std::abs(grad) < energyGradientCutoff_ )
							candidateCrossings[it].push_back(ibp);
					}
					if ( candidateCrossings[it].size() == 0 )
						throw std::runtime_error("No candidate for band continuation - cutoff too tight");
				}
				bool needOverlapCheck = false;
				for ( int it = 0;  it < targets.size(); ++it)
				{
					int iktir = kgrid.get_maps_red_to_irreducible()[targets[it]];
					needOverlapCheck = needOverlapCheck or (candidateCrossings[it].size() > 1);
					if ( candidateCrossings[it].size() == 1 )
					{
						bandOrder_[iktir*nB_ + ib] = candidateCrossings[it][0];
					}
				}
				if ( ! needOverlapCheck )
					continue;

				std::vector<std::vector<std::complex<float>>> reducWfcts;
				std::vector<int> npwPerK;
				wfct.generate_reducible_grid_wfcts(std::vector<int>{ibMap}, std::vector<int>{iksr}, reducWfcts, npwPerK);
				assert( npws == npwPerK[0] );
				std::fill(wfct1c.begin(), wfct1c.end(), std::complex<float>(0.0f));
				for (int ig = 0; ig < npws; ++ig)
					wfct1c[fms[ig*3+0]+fftMax[0]*(fms[ig*3+1]+fftMax[1]*fms[ig*3+2])] =
							std::conj(reducWfcts[0][ig]);


				for ( int it = 0;  it < targets.size(); ++it)
				{
					int iktr = targets[it];
					int iktir = kgrid.get_maps_red_to_irreducible()[targets[it]];
					if ( candidateCrossings[it].size() == 1 )
						continue;

					auto kt = kgrid.get_vector_direct(iktr);
					wfct.compute_Fourier_maps(kt, fftmaps);
					int npwt = fftmaps[0].size()/3;
					auto const & fmt = fftmaps[0];

					// here comes the part where we propagate the band convention from the k point iksr
					// to the target k point iktr
					std::map<float,int> ovlap;
					for ( int ibp : candidateCrossings[it])
					{
						assert((ibp >= 0) and (ibp < nB_));
						wfct.generate_reducible_grid_wfcts(std::vector<int>{ibp}, std::vector<int>{iktr}, reducWfcts, npwPerK);
						assert( npwt == npwPerK[0] );
						std::complex<float> overlap(0.0f);
						for (int ig = 0; ig < npwt; ++ig)
							overlap += wfct1c[fmt[ig*3+0]+fftMax[0]*(fmt[ig*3+1]+fftMax[1]*fmt[ig*3+2])]
											 *reducWfcts[0][ig];
						ovlap.insert(std::make_pair(std::abs(overlap), ibp));
					}
					auto c = ovlap.rbegin();
					bandOrder_[iktir*nB_ + ib] = c->second;
				}
			}
		}

//for ( auto l : linkedGrid )
//{
//	int ikr = l.first;
//	int ikir = kgrid.get_maps_red_to_irreducible()[ikr];
//	auto xyz = kgrid.get_reducible_to_xyz(ikr);
//	std::cout << xyz[0] << '\t'<< xyz[1] << '\t'<< xyz[2] << '\t'<< bandOrder_[ikir*nB_ + ib] << '\t' << std::endl;
//}
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
	auto const & d = kgrid.get_grid_dim();

	// represent both wavefunctions in the zone [0,1[
	std::vector<std::vector<int>> fftmaps;
	auto ks = kgrid.get_vector_direct(ikr1);
	auto kt = kgrid.get_vector_direct(ikr2);
	std::vector<int> g1zone(3);
	auto xyz1 = kgrid.get_reducible_to_xyz(ikr1);
	auto xyz2 = kgrid.get_reducible_to_xyz(ikr2);
	for ( int xi =0 ; xi < 3 ; ++xi )
	{
		double g = double(xyz1[xi])/double(d[xi]) - ks[xi];
		g1zone[xi] = std::floor( g + 0.5 );
		assert( std::abs(g - g1zone[xi]) < 1e-6 );
		g = double(xyz2[xi])/double(d[xi]) - kt[xi];
		assert( std::abs(g - std::floor( g + 0.5 )) < 1e-6 );
		g1zone[xi] = g1zone[xi] - std::floor( g + 0.5 );
	}
//
//	if ((g2zone[0] != g1zone[0])or((g2zone[1] != g1zone[1]))or(g2zone[2] != g1zone[2]))
//		std::cout << xyz1[0] <<'\t' << xyz1[1] <<'\t' << xyz1[2] <<
//		g1zone[0] <<'\t' << g1zone[1] <<'\t' << g1zone[2] <<'\t' <<
//		g2zone[0] <<'\t' << g2zone[1] <<'\t' << g2zone[2] <<'\t' << std::endl;


	wfct.compute_Fourier_maps(ks, fftmaps);
	int npws = fftmaps[0].size()/3;
	auto fms = fftmaps[0];

	auto fftMax = wfct.get_max_fft_dims();
	int npts = fftMax[0]*fftMax[1]*fftMax[2];
	std::vector<std::complex<float>> wfct1c(npts*bands1.size());

	std::vector<int> igocc;
	igocc.reserve(npws);
	std::vector<int> g(3);
	for (int ig = 0; ig < npws; ++ig)
	{
		for (int i=0 ;i <3 ; ++i)
			g[i] = fms[ig*3+i];
		Algorithms::FFTInterface::inplace_to_freq(g,fftMax);
		for (int i=0 ;i <3 ; ++i)
			g[i] -= g1zone[i];
		if ( (g[0] <= -fftMax[0]/2-fftMax[0]%2) or (g[0] > fftMax[0]/2)
				or (g[1] <= -fftMax[1]/2-fftMax[1]%2) or (g[1] > fftMax[1]/2)
				or (g[2] <= -fftMax[2]/2-fftMax[2]%2) or (g[2] > fftMax[2]/2))
			continue;
		Algorithms::FFTInterface::freq_to_inplace(g,fftMax);
		igocc.push_back(ig);
		for (int i=0 ;i <3 ; ++i)
			fms[ig*3+i] = g[i];
	}

	std::vector<std::vector<std::complex<float>>> reducWfcts;
	std::vector<int> npwPerK;
	wfct.generate_reducible_grid_wfcts(bands1, {ikr1}, reducWfcts, npwPerK);
	std::fill(wfct1c.begin(), wfct1c.end(), std::complex<float>(0.0f));
	for ( int ib1 : bands1)
		for (int ig : igocc)
		{
			wfct1c[fms[ig*3+0]+fftMax[0]*(fms[ig*3+1]+fftMax[1]*(fms[ig*3+2]+fftMax[2]*ib1))] =
					std::conj(reducWfcts[0][ig + npws*ib1]);
		}

	wfct.compute_Fourier_maps(kt, fftmaps);
	int npwt = fftmaps[0].size()/3;
	auto fmt = fftmaps[0];


	wfct.generate_reducible_grid_wfcts(bands2, {ikr2}, reducWfcts, npwPerK);

	overlaps.assign(nB_*nB_, 0.0f);
	for ( int ib1 = 0 ; ib1 < bands1.size(); ++ib1 )
		for (int ib2 = 0 ; ib2 < bands2.size(); ++ib2 )
		{
			if ( overlaps[bands1[ib1]+nB_*bands2[ib2]] != 0.0f )
				continue;

			std::complex<float> c(0.0f);
			for (int ig = 0 ; ig < npwt; ++ig)
				c += wfct1c[fmt[ig*3+0]+fftMax[0]*(fmt[ig*3+1]+fftMax[1]*(fmt[ig*3+2]+fftMax[2]*ib1))]
								 *reducWfcts[0][ig+npwt*ib2];
			overlaps[bands1[ib1]+nB_*bands2[ib2]] = std::abs(c);
			overlaps[bands2[ib2]+nB_*bands1[ib1]] = overlaps[bands1[ib1]+nB_*bands2[ib2]];
		}
}

int
BandOrderAnalysis::operator() (int ikir, int ib) const
{
	assert((ikir*nB_+ib >= 0) and (ikir*nB_+ib<bandOrder_.size()));
	return bandOrder_[ikir*nB_+ib];
}

} /* namespace ElectronicStructure */
} /* namespace elephon */
