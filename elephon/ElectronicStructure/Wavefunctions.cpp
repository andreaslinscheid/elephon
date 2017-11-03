/*	This file Wavefunctions.cpp is part of elephon.
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
 *  Created on: May 20, 2017
 *      Author: A. Linscheid
 */

#include "Wavefunctions.h"
#include "Algorithms/helperfunctions.hpp"
#include "Algorithms/FFTInterface.h"
#include <complex>
#include <iostream>
#include <set>

namespace elephon
{
namespace ElectronicStructure
{

void
Wavefunctions::initialize( LatticeStructure::RegularSymmetricGrid kgrid,
		std::shared_ptr<IOMethods::ElectronicStructureCodeInterface> wfctInterface)
{
	auto rootDir = wfctInterface->get_optns().get_root_dir();
	int numBands;
	wfctInterface->read_nBnd(rootDir, numBands);
	this->initialize(rootDir, std::move(kgrid), numBands, wfctInterface);
}

void
Wavefunctions::initialize(
		std::string rootDir,
		std::shared_ptr<IOMethods::ElectronicStructureCodeInterface> wfctInterface)
{
	LatticeStructure::RegularSymmetricGrid grid;
	wfctInterface->read_reciprocal_symmetric_grid(rootDir, grid);
	int numBands;
	wfctInterface->read_nBnd(rootDir, numBands);
	this->initialize(
			std::move(rootDir),
			std::move(grid),
			numBands,
			wfctInterface);
}

void
Wavefunctions::initialize(
		std::string rootDir,
		LatticeStructure::RegularSymmetricGrid kgrid,
		int numBands,
		std::shared_ptr<IOMethods::ElectronicStructureCodeInterface> wfctInterface)
{
	wfctInterface_ = wfctInterface;
	rootDir_ = std::move(rootDir);
	nBnd_ = numBands;
	tetraGrid_ = std::make_shared<LatticeStructure::TetrahedraGrid>();
	tetraGrid_->initialize(std::make_shared<LatticeStructure::RegularSymmetricGrid>(std::move(kgrid)));
}

void
Wavefunctions::generate_reducible_grid_wfcts(
		std::vector<int> const & bndIndices,
		std::vector<int> const & redKptIndices,
		std::vector< std::vector<std::complex<float>> > & wfcts,
		std::vector<int> & npwPerKAndSpin) const
{
	if ( not wfctInterface_ )
		throw std::invalid_argument("Calling uninitialized object");

	//First we make the connection reducible -> irreducible
	std::vector<int> irredKptIndices;
	auto grid = tetraGrid_->get_grid();
	grid->convert_reducible_irreducible(redKptIndices, irredKptIndices);

	//Find the unique irreducible indices and create a map that allows us to retrieve
	//the connection of input to these unique indices
	std::set<int> uniqIrrKptSet(irredKptIndices.begin(), irredKptIndices.end() );

	//Now we load the irreducible data
	std::vector<int> uniqIrrKpt(uniqIrrKptSet.begin(),	uniqIrrKptSet.end() );
	std::vector< std::vector< std::complex<float> > > irredWfcts;
	std::vector<int> npwPerKAndSpinIrred;
	wfctInterface_->read_wavefunctions(rootDir_,
			uniqIrrKpt,
			bndIndices,
			irredWfcts,
			npwPerKAndSpinIrred);

	//Allocate the storage for the reducible wave functions
	npwPerKAndSpin.resize(redKptIndices.size());
	for ( int ikred = 0; ikred < redKptIndices.size(); ++ikred )
	{
		auto it = uniqIrrKptSet.find( irredKptIndices[ikred] );
		int positionInUniqIrrKpts = std::distance(uniqIrrKptSet.begin(),it);
		assert( uniqIrrKpt[positionInUniqIrrKpts] == *it );
		npwPerKAndSpin[ikred] = npwPerKAndSpinIrred[ positionInUniqIrrKpts ];
	}

	wfcts.resize(redKptIndices.size());
	for ( int ikred = 0; ikred < redKptIndices.size(); ++ikred )
		wfcts[ikred].resize(npwPerKAndSpin[ikred]*bndIndices.size());

	//Make the connection of input indices and the position in the reducible grid
	std::map<int,int> reducibleToInputMap;
	for (int c = 0 ; c< redKptIndices.size(); ++c )
		reducibleToInputMap.insert( std::make_pair( redKptIndices[c], c ) );

	//Construct the Fourier mappings for all relevant irreducible kpoints
	std::vector<double> irredkpoints(3*uniqIrrKpt.size());
	for (int ik = 0 ; ik < uniqIrrKpt.size(); ++ik)
	{
		int ikIrred = uniqIrrKpt[ik];
		int isymId =  grid->get_maps_sym_irred_to_reducible()[ikIrred][ 0 ];
		assert( isymId == grid->get_symmetry().get_identity_index() );
		int ikRed = grid->get_maps_irreducible_to_reducible()[ikIrred][ isymId ];
		auto kIrred = grid->get_vector_direct(ikRed);
		std::copy(&kIrred[0],&kIrred[0]+3,&irredkpoints[ik*3]);
	}

	std::vector< std::vector<int> > irredFFTMap;
	wfctInterface_->compute_fourier_map( irredkpoints, irredFFTMap, grid->get_grid_prec());
	auto fftMaxDims = wfctInterface_->get_max_fft_dims();

	//And finally, we reconstruct the wavefunctions at reducible k points
	std::vector< std::vector<int> > irredToRed = grid->get_maps_irreducible_to_reducible();
	std::vector< std::vector<int> > symIrredToRed = grid->get_maps_sym_irred_to_reducible();
	for (int ik = 0 ; ik < uniqIrrKpt.size(); ++ik)
	{
		//ik enumerates the irreducible requested k points, ikir the corresponding irreducible index
		int ikir = uniqIrrKpt[ik];

		int npw = npwPerKAndSpinIrred[ik];

		//Construct a map that gives the consecutively ordered index
		// of the 3D G vector for a given plane wave index
		std::vector<int> pwToConseq(npw);
		assert( npw == irredFFTMap[ik].size()/3 );
		for ( int ipw = 0 ; ipw < npw; ++ipw)
		{
			int igx = irredFFTMap[ik][ipw*3+0];
			int igy = irredFFTMap[ik][ipw*3+1];
			int igz = irredFFTMap[ik][ipw*3+2];
			int conseq =  (igz*fftMaxDims[1]+igy)*fftMaxDims[0]+igx;
			pwToConseq[ipw] = conseq;
		}

		//using the irreducible vector, expand by symmetry - the requested vector must occur at least once
		for (int is = 0 ; is < static_cast<int>(symIrredToRed[ikir].size()); ++is)
		{
			//See if the reducible index appears in the requested set
			auto it = reducibleToInputMap.find( irredToRed[ikir][is] );
			if ( it == reducibleToInputMap.end() )
				continue;

			//construct a look up map to convert the rotated 3D vector index
			//back into a plane wave index
			int ikRed = irredToRed[ikir][ is ];
			auto kRed = grid->get_vector_direct(ikRed);
			std::vector< std::vector<int> > RedFFTMap;
			wfctInterface_->compute_fourier_map( kRed, RedFFTMap, grid->get_grid_prec());
			assert( npw == RedFFTMap[0].size()/3 );

			// Determine if there was a reciprocal lattice shift involved
			// This umklapp vector has to be added to the plane wave vector
			int isym = symIrredToRed[ikir][is];
			auto kIrred = grid->get_vector_direct(irredToRed[ikir][ grid->get_symmetry().get_identity_index() ]);
			auto ktransformNoPeriodic = kIrred;
			grid->get_symmetry().apply(isym, ktransformNoPeriodic, false);
			std::vector<int> umklappShift(3,0);
			for ( int i = 0 ; i < 3 ; ++i)
			{
				umklappShift[i] = std::floor(ktransformNoPeriodic[i] - kRed[i] + 0.5);
				assert( std::abs(umklappShift[i] - (ktransformNoPeriodic[i] - kRed[i])) < grid->get_grid_prec());
			}

			int inputLocation = it->second;
			if ( isym == grid->get_symmetry().get_identity_index() )
			{
				for ( int i = 0; i < int(bndIndices.size()); ++i)
					for ( int ipw = 0; ipw < npw; ++ipw)
						wfcts[inputLocation][ i*npw + ipw  ] =
								irredWfcts[ik][i*npw + ipw  ];
			}
			else//If this is not the identity symmetry we need to possibly
				// conjugate and/or multiply a phase and reorder the FFT mapping
			{
				bool inv = grid->get_symmetry().is_inversion( isym );
				bool symmorphic = grid->get_symmetry().is_symmorphic( isym );

				std::map<int,int> rotLookupMap;
				std::map<int,int>::const_iterator hint = rotLookupMap.end();

				for ( int ipw = 0 ; ipw < npw; ++ipw)
				{
					int igx = RedFFTMap[0][ipw*3+0];
					int igy = RedFFTMap[0][ipw*3+1];
					int igz = RedFFTMap[0][ipw*3+2];
					int conseq =  (igz*fftMaxDims[1]+igy)*fftMaxDims[0]+igx;
					hint = rotLookupMap.insert(hint, std::make_pair(conseq,ipw) );
				}

				//First we need to create a lookup table for the plane wave coefficients
				//associated to this wave function. The mappings are pre-computed and buffered
				this->fill_G_symmetry_buffer( isym );

				std::vector< std::complex<float> > phase;
				if ( not symmorphic )
					phase = std::vector< std::complex<float> >(npw);

				//For rotating the G vector, we need to 3D coordinates at this k point
				//Obtain the connection (index of pw) -> 3D G vector to rotate this wavefunction
				std::vector<int> G(3);
				std::vector<int> planeWaveLookup(npw,-1);
				for ( int ipw = 0 ; ipw < npw; ++ipw)
				{
					int cnsq = gSymBuffer_[isym][pwToConseq[ipw]];
					if ( (umklappShift[0] != 0) or (umklappShift[1] != 0) or (umklappShift[2] != 0) )
					{
						G[2] = cnsq/(fftMaxDims[1]*fftMaxDims[0]);
						G[1] = (cnsq-G[2]*(fftMaxDims[1]*fftMaxDims[0]))/fftMaxDims[0];
						G[0] = cnsq-(G[2]*fftMaxDims[1]+G[1])*fftMaxDims[0];
						Algorithms::FFTInterface::inplace_to_freq(G, fftMaxDims);
						for ( int i = 0 ; i < 3; ++i)
							G[i] += umklappShift[i];
						Algorithms::FFTInterface::freq_to_inplace(G, fftMaxDims);
						cnsq = (G[2]*fftMaxDims[1]+G[1])*fftMaxDims[0]+G[0];
					}

					auto ret = rotLookupMap.find(cnsq);
					if (ret == rotLookupMap.end())
					{
						G[2] = cnsq/(fftMaxDims[1]*fftMaxDims[0]);
						G[1] = (cnsq-G[2]*(fftMaxDims[1]*fftMaxDims[0]))/fftMaxDims[0];
						G[0] = cnsq-(G[2]*fftMaxDims[1]+G[1])*fftMaxDims[0];
						Algorithms::FFTInterface::inplace_to_freq(G, fftMaxDims);
						throw std::runtime_error(std::string("Could not locate rotated G vector ")
									+std::to_string(G[0]) + ", "+std::to_string(G[1]) + ", "+std::to_string(G[2]));
					}
					planeWaveLookup[ipw] = ret->second;

					assert( planeWaveLookup[ipw] < npw );
					if ( not symmorphic )
						phase[ipw] = phaseBuffer_[isym][pwToConseq[ipw]];
				}


				//If necessary compute the mapping
				//If the symmetry operation is non-symmorphic, we need to account
				// for the fractional translation by a phase shift in the plane wave coefficients
				if ( (not inv) and symmorphic )
					for ( int i = 0; i < int(bndIndices.size()); ++i)
						for ( int ipw = 0; ipw < npw; ++ipw)
							wfcts[inputLocation][ i*npw + planeWaveLookup[ipw]  ] =
									irredWfcts[ik][i*npw + ipw  ];

				if ( inv and symmorphic )
					for ( int i = 0; i < int(bndIndices.size()); ++i)
						for ( int ipw = 0; ipw < npw; ++ipw)
							wfcts[inputLocation][ i*npw + planeWaveLookup[ipw]  ] =
									std::conj(irredWfcts[ik][i*npw + ipw ]);

				if ( (not inv) and (not symmorphic) )
					for ( int i = 0; i < int(bndIndices.size()); ++i)
						for ( int ipw = 0; ipw < npw; ++ipw)
							wfcts[inputLocation][ i*npw + planeWaveLookup[ipw]  ] =
									phase[ipw]*irredWfcts[ik][i*npw + ipw ];

				if ( inv and (not symmorphic) )
					for ( int i = 0; i < int(bndIndices.size()); ++i)
						for ( int ipw = 0; ipw < npw; ++ipw)
							wfcts[inputLocation][ i*npw + planeWaveLookup[ipw]  ] =
									phase[ipw]*std::conj(irredWfcts[ik][i*npw + ipw ]);
			}
		}
	}
}

void
Wavefunctions::fill_G_symmetry_buffer(int isym) const
{
	auto symmetry = tetraGrid_->get_grid()->get_symmetry();
	auto fftDims = wfctInterface_->get_max_fft_dims();
	//Set the rotation buffer to the right size
	if ( int(gSymBuffer_.size()) != symmetry.get_num_symmetries() )
		gSymBuffer_ = std::vector<std::vector<int>>( symmetry.get_num_symmetries() );

	if ( int(gSymBuffer_[isym].size()) != fftDims[0]*fftDims[1]*fftDims[2] )
	{
		gSymBuffer_[isym].assign( fftDims[0]*fftDims[1]*fftDims[2] , -1);
		for ( int k = 0 ; k < fftDims[2]; ++k )
			for ( int j = 0 ; j < fftDims[1]; ++j )
				for ( int i = 0 ; i < fftDims[0]; ++i )
				{
					//By convention we store x as the fastest running index and z as the slowest
					int cnsq = (k*fftDims[1]+j)*fftDims[0]+i;
					std::vector<int> GRot = {i, j, k};
					Algorithms::FFTInterface::inplace_to_freq(GRot, fftDims);
					symmetry.rotate<int>(isym, GRot.begin(), GRot.end(), false);
					// We allow aliasing, because for uneven grids, rotated G vectors can exceed the fftDims.
					// These will not be referenced by plane waves anyhow.
					for ( int l = 0 ; l < 3 ; ++l)
					{
						GRot[l] = (GRot[l] <= -fftDims[l]/2-fftDims[l]%2 ? GRot[l] + fftDims[l] : GRot[l] );
						GRot[l] = (GRot[l] > fftDims[l]/2 ? GRot[l] - fftDims[l] : GRot[l] );
					}
					Algorithms::FFTInterface::freq_to_inplace(GRot, fftDims);
					int rotIndex = (GRot[2]*fftDims[1]+GRot[1])*fftDims[0]+GRot[0];
					assert( gSymBuffer_[isym][cnsq] < 0 );
					gSymBuffer_[isym][cnsq] = rotIndex;
				}
	}

	if ( not symmetry.is_symmorphic(isym) )
	{
		//Set the phase buffer to the right size
		if ( int(phaseBuffer_.size()) != symmetry.get_num_symmetries() )
			phaseBuffer_.resize( symmetry.get_num_symmetries() );

		if ( int(phaseBuffer_[isym].size()) != fftDims[0]*fftDims[1]*fftDims[2] )
		{
			std::vector<double> tau = symmetry.get_fractional_translation(isym);
			phaseBuffer_[isym].resize( fftDims[0]*fftDims[1]*fftDims[2] );
			for ( int k = 0 ; k < fftDims[2]; ++k )
				for ( int j = 0 ; j < fftDims[1]; ++j )
					for ( int i = 0 ; i < fftDims[0]; ++i )
					{
						std::vector<int> G = {i, j, k};
						Algorithms::FFTInterface::inplace_to_freq(G, fftDims);
						int consq = (k*fftDims[1]+j)*fftDims[0]+i;
						phaseBuffer_[isym][consq] = std::complex<float>(std::exp(
								 std::complex<double>(0,2.0*M_PI*(tau[0]*G[0]+tau[1]*G[1]+tau[2]*G[2])) ) );
					}
		}
	}
}

int
Wavefunctions::get_num_bands() const
{
	return nBnd_;
}

LatticeStructure::RegularSymmetricGrid const &
Wavefunctions::get_k_grid() const
{
	return *tetraGrid_->get_grid();
}

void
Wavefunctions::compute_Fourier_maps(std::vector<double> const & kvectors,
		std::vector<std::vector<int>> & fftMapsPerK) const
{
	wfctInterface_->compute_fourier_map( kvectors, fftMapsPerK, tetraGrid_->get_grid()->get_grid_prec());
}


void
Wavefunctions::generate_wfcts_at_arbitray_kp(
		std::vector<double> kList,
		std::vector<int> const & bandList,
		std::vector< std::vector< std::complex<float> > > & wfctsArbitrayKp,
		std::vector<std::vector<int>> & fftMapsArbitrayKp) const
{
	auto grid = tetraGrid_->get_grid();
	int nkA = kList.size()/3; //In the following A is for arbitrary location

	// Preparatory: Connect the input irregular grid points with the regular mesh.

	// Wave functions are on a regular grid and in G (reciprocal) space
	// First, find the tetrahedra where the k's are located and collect
	// all grid indices needed.
	std::vector<int> kAToCube;
	std::map<LatticeStructure::Tetrahedron, std::vector<int>> tetrasWithContainedKIndicesMap;
	for ( auto &kxi : kList )
		kxi -= std::floor(kxi+0.5);
	tetraGrid_->compute_grid_tetrahedra_surrounding_nongrid_points(
				kList,
				tetrasWithContainedKIndicesMap);
	// copy to a vector of pairs, since we need random access
	std::vector<std::pair<LatticeStructure::Tetrahedron, std::vector<int>>> tetrasWithContainedKIndices;
	tetrasWithContainedKIndices.reserve(tetrasWithContainedKIndicesMap.size());
	for (auto const & p : tetrasWithContainedKIndicesMap)
		tetrasWithContainedKIndices.push_back(p);

	// here: generate a list of all occurring reducible grid vectors in redIndices
	std::vector<int> npwPerKAllRedPts;
	std::vector<std::vector< std::complex<float> > > wfctAllRedPts;
	std::set<int> redIndicesSet;
	for ( auto const & t : tetrasWithContainedKIndices )
		redIndicesSet.insert(t.first.get_corner_indices().begin(), t.first.get_corner_indices().end());
	std::vector<int> redIndices(redIndicesSet.begin(),redIndicesSet.end());

	// load wavefunctions for this reducible indices.
	this->generate_reducible_grid_wfcts(
			bandList,
			redIndices,
			wfctAllRedPts,
			npwPerKAllRedPts);

	// Set up the FFT maps. Note this must be called after generate_reducible_grid_wfcts
	// since only then will the cutoffs and so on be read from the file.
	wfctsArbitrayKp.resize(nkA);
	this->compute_Fourier_maps(kList, fftMapsArbitrayKp);

	// Generate map from a corner point to the location in the array redIndices
	// and call it redIndicesTokptset
	std::vector<int> redIndicesTokptset( tetrasWithContainedKIndices.size()*4 );
	for ( int itetra = 0 ; itetra < tetrasWithContainedKIndices.size(); ++itetra )
		for ( int iCorner = 0 ; iCorner < 4; ++iCorner )
		{
			auto it = redIndicesSet.find( tetrasWithContainedKIndices[itetra].first.get_corner_indices()[iCorner] );
			assert( it != redIndicesSet.end() );
			redIndicesTokptset[itetra*4+iCorner] = std::distance(redIndicesSet.begin(),it);
			assert(redIndices[redIndicesTokptset[itetra*4+iCorner]] == *it);
		}
	// Preparatory done.

	// Step 1.: cast the wave functions at the corner points into a map for fast plane wave retrieval
	// 			Helper ft G Vector struct
	struct GVector {
		GVector( int x, int y, int z) : xi {x, y, z} { };
		int xi[3];
		bool operator<(GVector const & other) const
		{
			for ( int i = 0 ; i < 3; ++i)
				if ( xi[i] != other.xi[i] )
					return xi[i]<other.xi[i];
			return false;
		}
	};
	std::vector<double> redVectors(redIndices.size()*3);
	for ( int ikr = 0 ; ikr < redIndices.size() ; ++ikr )
	{
		auto v = grid->get_vector_direct( redIndices[ikr] );
		std::copy( v.begin(), v.end(), &redVectors[ikr*3] );
	}
	std::vector< std::vector<int> > redVecFFTMaps;
	this->compute_Fourier_maps(redVectors, redVecFFTMaps);

	auto fftMax = this->get_max_fft_dims();
	std::vector< std::map<GVector,int> >cmaps(redIndices.size());
	for ( int ikr = 0 ; ikr < redIndices.size() ; ++ikr )
		for ( int ipw = 0 ; ipw < redVecFFTMaps[ikr].size()/3; ++ipw)
		{
			int iGx = redVecFFTMaps[ikr][ipw*3+0];
			int iGy = redVecFFTMaps[ikr][ipw*3+1];
			int iGz = redVecFFTMaps[ikr][ipw*3+2];
			cmaps[ikr].insert( std::make_pair(GVector(iGx, iGy, iGz), ipw) );
		}

	std::vector<int> npwPerK(4);
	std::vector<std::vector< std::complex<float> > > wfctCornerPoints(4);
	for ( int itetra = 0 ; itetra < tetrasWithContainedKIndices.size(); ++itetra)
	{
	// Step 2.: Fetch the coefficients needed at this very k point.
	//			Coefficients that do not appear are set to zero.
		for ( int ikAC = 0 ; ikAC < tetrasWithContainedKIndices[itetra].second.size(); ++ikAC )
		{
			//ikA is the index of the k point in kList in this cube
			int ikA = tetrasWithContainedKIndices[itetra].second[ikAC];
			for ( int iCorner = 0 ; iCorner < 4 ; ++iCorner )
			{
				int ikr = redIndicesTokptset[itetra*4+iCorner];
				if ( fftMapsArbitrayKp[ikA] == redVecFFTMaps[ikr] )
				{
					wfctCornerPoints[iCorner] = wfctAllRedPts[ikr];
				}
				else
				{
					int npw = fftMapsArbitrayKp[ikA].size()/3;
					int npwC = redVecFFTMaps[ikr].size()/3;
					wfctCornerPoints[iCorner].resize( npw*bandList.size() );
					for ( int ib = 0 ; ib < bandList.size() ; ++ib )
					{
						for ( int ipw = 0 ; ipw < npw ; ++ipw )
						{
							int iGx = fftMapsArbitrayKp[ikA][ipw*3+0];
							int iGy = fftMapsArbitrayKp[ikA][ipw*3+1];
							int iGz = fftMapsArbitrayKp[ikA][ipw*3+2];
							auto it = cmaps[ikr].find(GVector(iGx, iGy, iGz));
							std::complex<float> val(0.0f);
							if ( it != cmaps[ikr].end() )
								val = wfctAllRedPts[ikr][ib*npwC+it->second];
							wfctCornerPoints[iCorner][ib*npw+ipw] = val;
						}
					}
				}
			}

	// Step 3.: Perform a linear interpolation. While the wavefunction need the actual 1BZ
	//			k vectors, here we have to map to the zone [0, 1[
			double kB[] = {	redVectors[redIndicesTokptset[itetra*4]*3 + 0],
							redVectors[redIndicesTokptset[itetra*4]*3 + 1],
							redVectors[redIndicesTokptset[itetra*4]*3 + 2]};
			std::vector<double> kv = {kList[ikA*3+0], kList[ikA*3+1], kList[ikA*3+2]};
			for ( int i = 0 ; i < 3 ; ++i )
			{
				kB[i] -= std::floor(kB[i]);
				kv[i] -= std::floor(kv[i]);
			}

	 		Algorithms::helperfunctions::interpolate_within_single_tetrahedron(
	 				kv,
					tetrasWithContainedKIndices[itetra].first,
					wfctCornerPoints,
					wfctsArbitrayKp[ikA]);
		}
	}
}

std::vector<int>
Wavefunctions::get_max_fft_dims() const
{
	return wfctInterface_->get_max_fft_dims();
}

} /* namespace ElectronicStructure */
} /* namespace elephon */
