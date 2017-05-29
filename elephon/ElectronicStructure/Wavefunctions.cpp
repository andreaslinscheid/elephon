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
#include <complex>
#include <iostream>

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
		std::vector< std::complex<float> > & wfcts,
		std::vector<int> & npwPerKAndSpin) const
{
	if ( not wfctInterface_ )
		throw std::invalid_argument("Calling uninitialized object");

	//First we make the connection reducible -> irreducible
	std::vector<int> irredKptIndices;
	grid_.convert_reducible_irreducible(redKptIndices,irredKptIndices);

	//Now we load the irreducible data
	std::set<int> uniqIrrKptSet(irredKptIndices.begin(),irredKptIndices.end());
	std::vector< std::complex<float> > irredWfcts;
	std::vector<int> npwPerKAndSpinIrred;
	wfctInterface_->read_wavefunctions(rootDir_,
			std::vector<int>(uniqIrrKptSet.begin(),	uniqIrrKptSet.end() ),
			bndIndices,irredWfcts,npwPerKAndSpinIrred);

	npwPerKAndSpin.resize(redKptIndices.size());
	for ( size_t ikred = 0; ikred < redKptIndices.size(); ++ikred )
		npwPerKAndSpin[ikred] = npwPerKAndSpinIrred[ irredKptIndices[ikred] ];

	//Create to lookups for reducible and irreducible wavefunction arrays.
	//irredPos[ik] will tell the beginning of the data for this k point in the
	//irreducible set of data values
	std::vector<int> irredPos( uniqIrrKptSet.size(), 0 );
	for ( int ikir = 1; ikir < int(npwPerKAndSpinIrred.size()); ++ikir )
		irredPos[ikir] = irredPos[ikir-1] + npwPerKAndSpinIrred[ikir-1]*bndIndices.size();

	std::vector<int> redPos( redKptIndices.size(), 0 );
	for ( int ikred = 1; ikred < int(redKptIndices.size()); ++ikred )
		redPos[ikred] = redPos[ikred-1] + npwPerKAndSpin[ikred-1]*bndIndices.size();

	//Allocate the storage for the reducible wave functions
	wfcts.resize( redPos.back()+npwPerKAndSpin.back()*bndIndices.size());

	//Make the connection of input indices and the position in the reducible grid
	std::map<int,int> reducibleToInputMap;
	for (int c = 0 ; c< int(redKptIndices.size()); ++c )
		reducibleToInputMap.insert( std::make_pair( redKptIndices[c], c ) );

	//Construct the Fourier mappings
	std::vector<double> irredkpoints(3*grid_.get_np_irred());
	for (int ikir = 0 ; ikir < grid_.get_np_irred(); ++ikir)
	{
		int ikRed = grid_.get_maps_irreducible_to_reducible()[ikir][ 0 ];
		auto kIrred = grid_.get_vector_direct(ikRed);
		std::copy(&kIrred[0],&kIrred[0]+3,&irredkpoints[ikir*3]);
	}
	std::vector< std::vector<int> > irredFFTMap;
	wfctInterface_->compute_fourier_map( irredkpoints, irredFFTMap);
	auto fftMaxDims = wfctInterface_->get_max_fft_dims();

	//And finally, we reconstruct the wavefunctions at reducible k points
	std::vector< std::vector<int> > irredToRed = grid_.get_maps_irreducible_to_reducible();
	std::vector< std::vector<int> > symIrredToRed = grid_.get_maps_sym_irred_to_reducible();
	for (int ikir = 0 ; ikir < grid_.get_np_irred(); ++ikir)
	{
		//See if this irreducible index appears in the requested set
		if ( uniqIrrKptSet.find( ikir ) == uniqIrrKptSet.end() )
			continue;

		int npw = npwPerKAndSpinIrred[ikir];

		//Construct a map that gives the consecutively ordered index
		// of the 3D G vector for a given plane wave index
		std::vector<int> pwToConseq(npwPerKAndSpinIrred[ikir]);
		assert( npw == int(irredFFTMap[ikir].size()/3) );
		for ( int ipw = 0 ; ipw < npw; ++ipw)
		{
			int igx = irredFFTMap[ikir][ipw*3+0];
			int igy = irredFFTMap[ikir][ipw*3+1];
			int igz = irredFFTMap[ikir][ipw*3+2];
			int conseq =  (igz*fftMaxDims[1]+igy)*fftMaxDims[0]+igx;
			pwToConseq[ipw] = conseq;
		}

		for (int is = 0 ; is < static_cast<int>(symIrredToRed[ikir].size()); ++is)
		{
			//See if the reducible index appears in the requested set
			auto it = reducibleToInputMap.find( irredToRed[ikir][is] );
			if ( it == reducibleToInputMap.end() )
				continue;

			//construct a look up map to convert the rotated 3D vector index
			//back into a plane wave index
			int ikRed = grid_.get_maps_irreducible_to_reducible()[ikir][ is ];
			auto kRed = grid_.get_vector_direct(ikRed);
			std::vector< std::vector<int> > RedFFTMap;
			wfctInterface_->compute_fourier_map( kRed, RedFFTMap);
			assert( npw == int(RedFFTMap[0].size()/3) );

			int inputLocation = it->second;

			int isym = symIrredToRed[ikir][is];
			if ( isym == grid_.get_symmetry().get_identity_index() )
			{
				for ( int i = 0; i < int(bndIndices.size()); ++i)
					for ( int ipw = 0; ipw < npw; ++ipw)
						wfcts[ redPos[inputLocation] + i*npw + ipw  ] =
								irredWfcts[ irredPos[ikir] + i*npw + ipw  ];
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
					auto ret = rotLookupMap.find(gSymBuffer_[isym][pwToConseq[ipw]]);
					if ( ret == rotLookupMap.end() )
						throw std::runtime_error( " Could not located rotated G vector " );

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
							wfcts[ redPos[inputLocation] + i*npw + planeWaveLookup[ipw]  ] =
									irredWfcts[ irredPos[ikir] + i*npw + ipw  ];

				if ( inv and symmorphic )
					for ( int i = 0; i < int(bndIndices.size()); ++i)
						for ( int ipw = 0; ipw < npw; ++ipw)
							wfcts[ redPos[inputLocation] + i*npw + planeWaveLookup[ipw]  ] =
									std::conj(irredWfcts[ irredPos[ikir] + i*npw + ipw ]);

				if ( (not inv) and (not symmorphic) )
					for ( int i = 0; i < int(bndIndices.size()); ++i)
						for ( int ipw = 0; ipw < npw; ++ipw)
							wfcts[ redPos[inputLocation] + i*npw + planeWaveLookup[ipw]  ] =
									phase[ipw]*irredWfcts[ irredPos[ikir] + i*npw + ipw ];

				if ( inv and (not symmorphic) )
					for ( int i = 0; i < int(bndIndices.size()); ++i)
						for ( int ipw = 0; ipw < npw; ++ipw)
							wfcts[ redPos[inputLocation] + i*npw + planeWaveLookup[ipw]  ] =
									phase[ipw]*std::conj(irredWfcts[ irredPos[ikir] + i*npw + ipw ]);
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

	if ( int(gSymBuffer_[isym].size()) != fftDims[0]*fftDims[1]*fftDims[2] )
	{
		gSymBuffer_[isym] = std::vector<int>( fftDims[0]*fftDims[1]*fftDims[2] );
		auto S = grid_.get_symmetry().get_sym_op(isym);
		for ( int k = 0 ; k < fftDims[2]; ++k )
			for ( int j = 0 ; j < fftDims[1]; ++j )
				for ( int i = 0 ; i < fftDims[0]; ++i )
				{
					//By convention we store x as the fastest running index and z as the slowest
					int indexOrig = (k*fftDims[1]+j)*fftDims[0]+i;
					int gx = i < fftDims[0]/2 ? i : i - fftDims[0];
					int gy = j < fftDims[1]/2 ? j : j - fftDims[1];
					int gz = k < fftDims[2]/2 ? k : k - fftDims[2];
					int GRot[3] = {	S.ptgroup[0*3+0]*gx+S.ptgroup[0*3+1]*gy+S.ptgroup[0*3+2]*gz,
									S.ptgroup[1*3+0]*gx+S.ptgroup[1*3+1]*gy+S.ptgroup[1*3+2]*gz,
									S.ptgroup[2*3+0]*gx+S.ptgroup[2*3+1]*gy+S.ptgroup[2*3+2]*gz};
					int ri = GRot[0] < 0 ? GRot[0] + fftDims[0]: GRot[0];
					assert( (ri >= 0) && (ri < fftDims[0]) );
					int rj = GRot[1] < 0 ? GRot[1] + fftDims[1]: GRot[1];
					assert( (rj >= 0) && (rj < fftDims[1]) );
					int rk = GRot[2] < 0 ? GRot[2] + fftDims[2]: GRot[2];
					assert( (rk >= 0) && (rk < fftDims[2]) );
					int rotIndex = (rk*fftDims[1]+rj)*fftDims[0]+ri;
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

} /* namespace ElectronicStructure */
} /* namespace elephon */
