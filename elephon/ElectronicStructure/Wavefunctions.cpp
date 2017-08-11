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
#include "Algorithms/TrilinearInterpolation.h"
#include <complex>
#include <iostream>
#include <set>

namespace elephon
{
namespace ElectronicStructure
{

void
Wavefunctions::initialize( double gridPrec,
		std::string rootDir,
		std::shared_ptr<IOMethods::ElectronicStructureCodeInterface> wfctInterface)
{
	wfctInterface_ = std::move( wfctInterface );
	rootDir_ = std::move(rootDir);

	//Read and set up the symmetry module and the lattice module
	LatticeStructure::LatticeModule lattice;
	LatticeStructure::Symmetry sym;
	std::vector<LatticeStructure::Atom> atoms;
	wfctInterface_->read_cell_paramters(rootDir_, gridPrec, grid_, lattice, atoms, sym);
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
	grid_.convert_reducible_irreducible(redKptIndices, irredKptIndices);

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
		int isymId =  grid_.get_maps_sym_irred_to_reducible()[ikIrred][ 0 ];
		assert( isymId == grid_.get_symmetry().get_identity_index() );
		int ikRed = grid_.get_maps_irreducible_to_reducible()[ikIrred][ isymId ];
		auto kIrred = grid_.get_vector_direct(ikRed);
		std::copy(&kIrred[0],&kIrred[0]+3,&irredkpoints[ik*3]);
	}

	std::vector< std::vector<int> > irredFFTMap;
	wfctInterface_->compute_fourier_map( irredkpoints, irredFFTMap, grid_.get_grid_prec());
	auto fftMaxDims = wfctInterface_->get_max_fft_dims();

	//And finally, we reconstruct the wavefunctions at reducible k points
	std::vector< std::vector<int> > irredToRed = grid_.get_maps_irreducible_to_reducible();
	std::vector< std::vector<int> > symIrredToRed = grid_.get_maps_sym_irred_to_reducible();
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
			auto kRed = grid_.get_vector_direct(ikRed);
			std::vector< std::vector<int> > RedFFTMap;
			wfctInterface_->compute_fourier_map( kRed, RedFFTMap, grid_.get_grid_prec());
			assert( npw == RedFFTMap[0].size()/3 );

			//Determine if there was a reciprocal lattice shift involved
			//This umklapp vector has to be added to the plane wave vector
			int isym = symIrredToRed[ikir][is];
			auto kIrred = grid_.get_vector_direct(irredToRed[ikir][ 0 ]);
			auto ktransformNoPeriodic = kIrred;
			grid_.get_symmetry().apply(isym, ktransformNoPeriodic, false);
			std::vector<int> umklappShift(3,0);
			for ( int i = 0 ; i < 3 ; ++i)
			{
				umklappShift[i] = std::floor(ktransformNoPeriodic[i] - kRed[i] + 0.5);
				assert( std::abs(umklappShift[i] - (ktransformNoPeriodic[i] - kRed[i])) < grid_.get_grid_prec());
			}

			int inputLocation = it->second;
			if ( isym == grid_.get_symmetry().get_identity_index() )
			{
				for ( int i = 0; i < int(bndIndices.size()); ++i)
					for ( int ipw = 0; ipw < npw; ++ipw)
						wfcts[inputLocation][ i*npw + ipw  ] =
								irredWfcts[ik][i*npw + ipw  ];
			}
			else//If this is not the identity symmetry we need to possibly
				// conjugate and/or multiply a phase and reorder the FFT mapping
			{
				bool inv = grid_.get_symmetry().is_inversion( isym );
				bool symmorphic = grid_.get_symmetry().is_symmorphic( isym );

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
				std::vector<int> planeWaveLookup(npw,-1);
				for ( int ipw = 0 ; ipw < npw; ++ipw)
				{
					int cnsq = gSymBuffer_[isym][pwToConseq[ipw]];
					if ( (umklappShift[0] != 0) or (umklappShift[1] != 0) or (umklappShift[2] != 0) )
					{
						int iGz = cnsq/(fftMaxDims[1]*fftMaxDims[0]);
						int iGy = (cnsq-iGz*(fftMaxDims[1]*fftMaxDims[0]))/fftMaxDims[0];
						int iGx = cnsq-(iGz*fftMaxDims[1]+iGy)*fftMaxDims[0];
						iGx = iGx < fftMaxDims[0]/2 ?  iGx : iGx - fftMaxDims[0];
						iGy = iGy < fftMaxDims[1]/2 ?  iGy : iGy - fftMaxDims[1];
						iGz = iGz < fftMaxDims[2]/2 ?  iGz : iGz - fftMaxDims[2];
						iGx -= umklappShift[0];
						iGy -= umklappShift[1];
						iGz -= umklappShift[2];
						iGx = iGx < 0 ?  iGx + fftMaxDims[0] : iGx;
						iGy = iGy < 0 ?  iGy + fftMaxDims[1] : iGy;
						iGz = iGz < 0 ?  iGz + fftMaxDims[2] : iGz;
						cnsq = (iGz*fftMaxDims[1]+iGy)*fftMaxDims[0]+iGx;
					}

					auto ret = rotLookupMap.find(cnsq);
					if (ret == rotLookupMap.end())
					{
						int iGz = cnsq/(fftMaxDims[1]*fftMaxDims[0]);
						int iGy = (cnsq-iGz*(fftMaxDims[1]*fftMaxDims[0]))/fftMaxDims[0];
						int iGx = cnsq-(iGz*fftMaxDims[1]+iGy)*fftMaxDims[0];
						iGx = iGx < fftMaxDims[0]/2 ?  iGx : iGx - fftMaxDims[0];
						iGy = iGy < fftMaxDims[1]/2 ?  iGy : iGy - fftMaxDims[1];
						iGz = iGz < fftMaxDims[2]/2 ?  iGz : iGz - fftMaxDims[2];
						throw std::runtime_error(std::string("Could not located rotated G vector ")
									+std::to_string(iGx) + ", "+std::to_string(iGy) + ", "+std::to_string(iGz));
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
	auto fftDims = wfctInterface_->get_max_fft_dims();
	//Set the rotation buffer to the right size
	if ( int(gSymBuffer_.size()) != grid_.get_symmetry().get_num_symmetries() )
		gSymBuffer_ = std::vector<std::vector<int>>( grid_.get_symmetry().get_num_symmetries() );

	std::vector<int> GRot(3);
	if ( int(gSymBuffer_[isym].size()) != fftDims[0]*fftDims[1]*fftDims[2] )
	{
		gSymBuffer_[isym] = std::vector<int>( fftDims[0]*fftDims[1]*fftDims[2] , -1);
		auto S = grid_.get_symmetry().get_sym_op(isym);
		for ( int k = 0 ; k < fftDims[2]; ++k )
			for ( int j = 0 ; j < fftDims[1]; ++j )
				for ( int i = 0 ; i < fftDims[0]; ++i )
				{
					//By convention we store x as the fastest running index and z as the slowest
					int indexOrig = (k*fftDims[1]+j)*fftDims[0]+i;
					GRot[0] = i < fftDims[0]/2 ? i : i - fftDims[0];
					GRot[1] = j < fftDims[1]/2 ? j : j - fftDims[1];
					GRot[2] = k < fftDims[2]/2 ? k : k - fftDims[2];
					S.rotate(GRot);
					int ri = (GRot[0] < 0 ? GRot[0] + fftDims[0]: GRot[0]);
					int rj = (GRot[1] < 0 ? GRot[1] + fftDims[1]: GRot[1]);
					int rk = (GRot[2] < 0 ? GRot[2] + fftDims[2]: GRot[2]);
					ri = (ri == fftDims[0] ? 0 : ri);
					rj = (rj == fftDims[1] ? 0 : rj);
					rk = (rk == fftDims[2] ? 0 : rk);
					assert( (ri >= 0) && (ri < fftDims[0]) );
					assert( (rj >= 0) && (rj < fftDims[1]) );
					assert( (rk >= 0) && (rk < fftDims[2]) );
					int rotIndex = (rk*fftDims[1]+rj)*fftDims[0]+ri;
					assert( gSymBuffer_[isym][indexOrig] < 0 );
					gSymBuffer_[isym][indexOrig] = rotIndex;
				}
	}

	if ( not grid_.get_symmetry().is_symmorphic(isym) )
	{
		//Set the phase buffer to the right size
		if ( int(phaseBuffer_.size()) != grid_.get_symmetry().get_num_symmetries() )
			phaseBuffer_.resize( grid_.get_symmetry().get_num_symmetries() );

		if ( int(phaseBuffer_[isym].size()) != fftDims[0]*fftDims[1]*fftDims[2] )
		{
			std::vector<double> tau = grid_.get_symmetry().get_fractional_translation(isym);
			phaseBuffer_[isym].resize( fftDims[0]*fftDims[1]*fftDims[2] );
			for ( int k = 0 ; k < fftDims[2]; ++k )
				for ( int j = 0 ; j < fftDims[1]; ++j )
					for ( int i = 0 ; i < fftDims[0]; ++i )
					{
						int gx = i < fftDims[0]/2 ? i : i - fftDims[0];
						int gy = j < fftDims[1]/2 ? j : j - fftDims[1];
						int gz = k < fftDims[2]/2 ? k : k - fftDims[2];
						int consq = (k*fftDims[1]+j)*fftDims[0]+i;
						phaseBuffer_[isym][consq] = std::complex<float>(std::exp(
								 std::complex<double>(0,2.0*M_PI*(tau[0]*gx+tau[1]*gy+tau[2]*gz)) ) );
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
	return grid_;
}

void
Wavefunctions::compute_Fourier_maps(std::vector<double> const & kvectors,
		std::vector<std::vector<int>> & fftMapsPerK) const
{
	wfctInterface_->compute_fourier_map( kvectors, fftMapsPerK, grid_.get_grid_prec());
}


void
Wavefunctions::generate_wfcts_at_arbitray_kp(
		std::vector<double> kList,
		std::vector<int> const & bandList,
		std::vector< std::vector< std::complex<float> > > & wfctsArbitrayKp,
		std::vector<std::vector<int>> & fftMapsArbitrayKp) const
{
	int nkA = kList.size()/3; //In the following A is for arbitrary location

	//reserve space for the results
	wfctsArbitrayKp.resize(nkA);
	fftMapsArbitrayKp.resize(nkA);

	//Wave functions are on a regular grid and in G (reciprocal) space
	//First, load the cubes of wave functions for linear interpolation
	std::vector<int> kAToCube;
	std::vector<LatticeStructure::RegularSymmetricGrid::GridCube> gridCubes;
	for ( auto &kxi : kList )
		kxi -= std::floor(kxi);
	grid_.compute_grid_cubes_surrounding_nongrid_points(
				kList, kAToCube, gridCubes);

	//Generate a list of all occurring reducible grid vectors.
	std::vector<int> npwPerKAllRedPts;
	std::vector<std::vector< std::complex<float> > > wfctAllRedPts;
	std::set<int> redIndicesSet;
	for ( auto c : gridCubes )
		redIndicesSet.insert( c.cornerIndices_.begin(),  c.cornerIndices_.end() );
	std::vector<int> redIndices(redIndicesSet.begin(),redIndicesSet.end());
	this->generate_reducible_grid_wfcts(
			bandList,
			redIndices,
			wfctAllRedPts,
			npwPerKAllRedPts);

	//Generate map from a corner point to the location in the array redIndices
	std::vector<int> redIndicesTokptset( gridCubes.size()*8 );
	for ( int ic = 0 ; ic < gridCubes.size(); ++ic )
		for ( int iCorner = 0 ; iCorner < 8; ++iCorner )
		{
			auto it = redIndicesSet.find( gridCubes[ic].cornerIndices_[iCorner] );
			assert( it != redIndicesSet.end() );
			redIndicesTokptset[ic*8+iCorner] = std::distance(redIndicesSet.begin(),it);
			assert(redIndices[redIndicesTokptset[ic*8+iCorner]] == *it);
		}

	Algorithms::TrilinearInterpolation interpol( grid_.view_bare_grid() );

	std::vector<int> npwPerK(8);
	std::vector<std::vector< std::complex<float> > > wfctCornerPoints(8);
	int iIrrKCounter = 0;
	for ( int ic = 0 ; ic < gridCubes.size(); ++ic)
	{
		//fetch the wave functions at the corner points
		for ( int iCorner = 0 ; iCorner < 8; ++iCorner )
		{
			npwPerK[iCorner] = npwPerKAllRedPts[redIndicesTokptset[ic*8+iCorner]];
			wfctCornerPoints[iCorner] = wfctAllRedPts[redIndicesTokptset[ic*8+iCorner]];
		}

		//Possibly adjust the size to the largest common set of plane waves
		std::vector<int> commonFFTMap;
		//First get the fft maps for the corner points
		std::vector<double> cornerVectors(8*3);
		for ( int iCorner = 0 ; iCorner < 8 ; ++iCorner )
		{
			auto v = grid_.get_vector_direct( gridCubes[ic].cornerIndices_[iCorner] );
			std::copy( v.begin(), v.end(), &cornerVectors[iCorner*3] );
		}
		std::vector< std::vector<int> > cornerVecFFTMaps;
		this->compute_Fourier_maps(cornerVectors, cornerVecFFTMaps);

		//If they do not all coincide, we must generate a new FFT map the includes all occurring plane waves
		bool newMap = false;
		for ( int iCorner = 1 ; iCorner < 8 ; ++iCorner )
			if ( cornerVecFFTMaps[0] != cornerVecFFTMaps[iCorner] )
			{
				newMap = true;
				break;
			}

		if ( newMap )
		{
			//Helper ft G Vector struct
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
			//get the respective min/max values among the corner points in x,y and z
			auto fftMax = this->get_max_fft_dims();
			std::vector<int> max{0,0,0};
			auto min = fftMax;
			std::vector< std::map<GVector,int> >cmaps(8);
			for ( int iCorner = 0 ; iCorner < 8 ; ++iCorner )
				for ( int ipw = 0 ; ipw < cornerVecFFTMaps[iCorner].size(); ++ipw)
				{
					int iGx = cornerVecFFTMaps[iCorner][ipw*3+0];
					int iGy = cornerVecFFTMaps[iCorner][ipw*3+1];
					int iGz = cornerVecFFTMaps[iCorner][ipw*3+2];
					iGx = iGx < fftMax[0]/2 ? iGx : iGx - fftMax[0];
					iGy = iGy < fftMax[1]/2 ? iGy : iGy - fftMax[1];
					iGz = iGz < fftMax[2]/2 ? iGz : iGz - fftMax[2];
					max[0] = std::max(max[0],iGx);
					max[1] = std::max(max[1],iGy);
					max[2] = std::max(max[2],iGz);
					min[0] = std::min(min[0],iGx);
					min[1] = std::min(min[1],iGy);
					min[2] = std::min(min[2],iGz);
					cmaps[iCorner].insert( std::make_pair(GVector(iGx,iGy,iGz),ipw) );
				}

			commonFFTMap.reserve(3*fftMax[0]*fftMax[1]*fftMax[2]);

			std::vector< std::vector<std::complex<float> > > wfctCornerPointsBuffer(8);
			for ( int iCorner = 0 ; iCorner < 8 ; ++iCorner )
				wfctCornerPointsBuffer[iCorner].reserve(fftMax[0]*fftMax[1]*fftMax[2]);

			//Now expand each corner wavefunction with zeros such that all corners are
			//given on the same set of G vector - even though some coefficients are zero of cause
			for ( int iGz = min[2] ; iGz <= max[2]; ++iGz)
				for ( int iGy = min[1] ; iGy <= max[1]; ++iGy)
					for ( int iGx = min[0] ; iGx <= max[0]; ++iGx)
					{
						std::vector<int> occursAtPWIndex(8,-1);
						for ( int iCorner = 0 ; iCorner < 8 ; ++iCorner )
						{
							auto ret = cmaps[iCorner].find( GVector(iGx,iGy,iGz) );
							if ( ret != cmaps[iCorner].end() )
								occursAtPWIndex[iCorner] = ret->second ;
						}

						if ( *std::max_element(occursAtPWIndex.begin(),occursAtPWIndex.end()) == -1 )
							//This plane wave does not occur anywhere ...
							continue;

						// Since this plane wave occurs at some of the corner wavefunctions it goes into
						// the map
						int iGxf = iGx < 0 ? iGx + fftMax[0] : iGx;
						int iGyf = iGy < 0 ? iGy + fftMax[1] : iGy;
						int iGzf = iGz < 0 ? iGz + fftMax[2] : iGz;
						assert((iGxf >= 0) && (iGxf < fftMax[0]) && (iGyf >= 0) && (iGyf < fftMax[1])
							&& (iGzf >= 0) && (iGzf < fftMax[2]));
						commonFFTMap.push_back(iGxf);
						commonFFTMap.push_back(iGyf);
						commonFFTMap.push_back(iGzf);

						//Elements that do not occur are set to zero.
						for ( int iCorner = 0 ; iCorner < 8 ; ++iCorner )
							for ( int ib = 0 ; ib < bandList.size() ; ++ib)
							{
								int indexOldPWSet = occursAtPWIndex[iCorner]*bandList.size()+ib;
								wfctCornerPointsBuffer[iCorner].push_back(
										occursAtPWIndex[iCorner] < 0 ? std::complex<float>(0) :
												wfctCornerPoints[iCorner][indexOldPWSet] );
							}
					}

			for ( int iCorner = 0 ; iCorner < 8 ; ++iCorner )
				wfctCornerPoints[iCorner] = wfctCornerPointsBuffer[iCorner];
		}
		else
			commonFFTMap = cornerVecFFTMaps[0];

		int npw = wfctCornerPoints[0].size()/bandList.size();
		int nB = bandList.size();

		//v0 is the vector pointing to the lowest x,y and z corner point in direct coordinates
		auto v0 = grid_.get_vector_direct( gridCubes[ic].cornerIndices_[0] );
		std::vector<double> relIrregularKCoords( gridCubes[ic].containedIrregularPts_.size()*3 );
		for ( int ikAC = 0 ; ikAC < gridCubes[ic].containedIrregularPts_.size(); ++ikAC )
		{
			int ikA = gridCubes[ic].containedIrregularPts_[ikAC];
			std::copy(&kList[ ikA*3 ], &kList[ ikA*3 ]+3 , &relIrregularKCoords[ikAC*3] );
		}

		std::vector< std::complex<float> > buffer;
 		interpol.interpolate_within_single_cube(
				relIrregularKCoords,
				wfctCornerPoints, buffer);

 		for ( int ikAt = 0 ; ikAt < gridCubes[ic].containedIrregularPts_.size(); ++ikAt)
 		{
 			int i = iIrrKCounter+ikAt; //the total index in the list of irregular k points
 	 		wfctsArbitrayKp[i].insert( wfctsArbitrayKp[i].end(), &buffer[ikAt*npw*nB] , &buffer[ikAt*npw*nB] + npw*nB );
 	 		fftMapsArbitrayKp[i] = commonFFTMap;
 		}
		iIrrKCounter += gridCubes[ic].containedIrregularPts_.size();
	}
}

std::vector<int>
Wavefunctions::get_max_fft_dims() const
{
	return wfctInterface_->get_max_fft_dims();
}

} /* namespace ElectronicStructure */
} /* namespace elephon */
