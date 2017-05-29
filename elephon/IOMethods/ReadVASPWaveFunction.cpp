/*	This file ReadVASPWaveFunction.cpp is part of elephon.
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
 *  Created on: Apr 24, 2017
 *      Author: A. Linscheid
 */

#include <IOMethods/ReadVASPWaveFunction.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <assert.h>

namespace elephon
{
namespace IOMethods
{

ReadVASPWaveFunction::ReadVASPWaveFunction()
{
	// TODO Auto-generated constructor stub

}

ReadVASPWaveFunction::~ReadVASPWaveFunction()
{
	// TODO Auto-generated constructor stub

}

void ReadVASPWaveFunction::prepare_wavecar(std::string filename)
{
	if ( filename_.compare( filename ) == 0 )
		return;

	filename_ = std::move( filename );

	if ( wavecarfile_.is_open() )
		wavecarfile_.close();

	//General comments from looking at the VASP source
	//
	// 1) All integer values are stored as floating point numbers
	//
	// 2) The precision VASP uses in this file is
	// 		selected_real_kind( 10 ) which means that it is the
	//		smallest floating point number capable of storing 10 digits
	//		the processor has to offer. This is double on amd64.
	//		typedef double VASPDprec;
	//		typedef float VASPSprec;
	//TODO Make sure this works on other platforms too.
	//
	// 3) The make file for intel seems to inlcude '-assume byterecl'
	//		thus we assume at this point that a record block is one byte.
	//		The default for intel is 4 bytes. We don't have a reliable method to determine
	//		what block size actually was and if somebody compiles VASP without
	//		'-assume byterecl' it will break this code.
	//
	// 4) In more modern versions, VASP can put the data over a given record.
	//		This means the file is not segmented any more and we have lost
	//		the only advantage Fortran 'direct' access really gives us.
	//		The flag is wavefuncSpanRecords_ and is inferred from the VASP version identifier
	//		This means in this case we need to walk through the entire file first and
	//		save where the data for (ik,ibnd,ispin) begins.
	//
	//		Fortran I/O is *TRUELY* useless.

	wavecarfile_.open( filename_.c_str(), std::ios::in | std::ios::binary);
	if ( ! wavecarfile_.good() )
		throw std::runtime_error( std::string("File not readable: " + filename_));
	wavecarfile_.seekg(0, std::ios::end);
	total_size_ = wavecarfile_.tellg();
	wavecarfile_.seekg(0);
	//First, we read the header which is always small.
	//The first three numbers are
	// Rrecl - the Fortran record length
	// Rnspin - the number of spins
	// Rversion - an integer signaling the version of VASP that wrote this file

	auto nint = [](VASPDprec a) { return static_cast<int>(std::floor(a+0.5)); };

	//First we check the actual record length, reset our buffer and the view in terms of floating point values.
	size_t saveSize = (1024 < total_size_ ? 1024: total_size_); // should be enough to read the first line header
	std::vector<char> buffer(saveSize);
	wavecarfile_.read( &buffer[0] , saveSize);
	VASPDprec const * viewAsFloat = reinterpret_cast<VASPDprec const *>( &buffer[0] );
	recl_ = static_cast<std::size_t>( std::floor(viewAsFloat[0]+0.5) );
	nspin_ = nint(viewAsFloat[1]);

	int version = nint(viewAsFloat[2]);
	if ( version == 53300 ) // "VASP.5.3 WAVECAR encountered"
	{
		wavefuncDouble_  = false;
		spanRecords_ = true;
	}
	else if( version == 53310) //"VASP.5.3 double precision WAVECAR encountered"
	{
		wavefuncDouble_ = true;
		spanRecords_ = true;
	}
	else if ( version == 45200 ) //"VASP 5.4.1 WAVECAR encountered"
	{
		wavefuncDouble_ = false;
		spanRecords_ = false;
	}
	else //"double precision WAVECAR encountered"
	{
		wavefuncDouble_ = true;
		spanRecords_ = false;
	}

	//go to the second record locator in the file
	//This contains the #kpts, #bands and then the 9 double values for the lattice matrix
	//NOTE: This matrix is defined in the same order we use in elephon
	wavecarfile_.seekg(recl_);
	wavecarfile_.read( &buffer[0], saveSize);

	int nkptsVASP =nint(viewAsFloat[0]);
	nBndsVASP_=nint(viewAsFloat[1]);
	ecutoff_ = viewAsFloat[2];

	std::vector<double> latticeMatrix(&(viewAsFloat[3]),&(viewAsFloat[3])+9);
	lattice_.initialize( latticeMatrix );

	this->set_up_fourier_max( );

	//Now we walk through the file and fetch the band energies and k points.
	kpoints_ = std::vector<double>( nkptsVASP*3 );
	energies_ = std::vector<double>( nBndsVASP_*nkptsVASP );
	kptSpinPosToFile_ = std::vector<std::size_t>( nspin_* nkptsVASP );
	npwSpinKpt_ = std::vector<int>(nkptsVASP*nspin_);

	//Reset the buffer (and its view) to the size we need to read at a time
	std::size_t byteNumBandsHeader = static_cast<std::size_t>(4+nBndsVASP_*3)*sizeof(VASPDprec);
	std::size_t numBytesRead = std::max(recl_,byteNumBandsHeader);
	buffer = std::vector<char>(numBytesRead);
	viewAsFloat = reinterpret_cast<VASPDprec const *>( &buffer[0] );

	std::size_t irecord = 2;
	for (int ispin = 0 ; ispin < nspin_; ++ispin)
		for (int ik = 0 ; ik < nkptsVASP; ++ik)
		{
			 // We only care about the first values with the band energies
			assert( (irecord +1)* recl_ <= total_size_ );
			wavecarfile_.seekg( irecord * recl_ );
			wavecarfile_.read( &buffer[0], numBytesRead  );

			//Apart from the normal increase in irecord above, we must add addition advances due to
			//the possibility that the data may span records.
			irecord += (1 + this->num_records_spanned( numBytesRead ) );
			//Now comes the wave function data.
			//Keep track of where we are in the file for a future wave function reading
			kptSpinPosToFile_[ik*nspin_+ispin] = irecord;

			//Position 1 is the number of plane wave
			npwSpinKpt_[ik*nspin_+ispin] = nint(viewAsFloat[0]);
			//Position 2-4 is the k vector.
			std::copy(&viewAsFloat[1],&viewAsFloat[1]+3,&kpoints_[ik*3]);

			size_t startBandsLoop = 4;
			//now we loop the bands and copy the energies
			for (int ibnd = 0 ; ibnd < nBndsVASP_; ++ibnd)
			{
				//First is the real part of the energies
				energies_[(ik*nBndsVASP_+ibnd)*nspin_+ispin] = viewAsFloat[startBandsLoop];
				//Then comes the complex part ... for whatever reason
				//Now comes the occupation : makes 3 double values per band
				startBandsLoop += 3;
			}

			//In case spanRecords is false, the wave function to this spin/band/kpt are
			//all in this same record and we are done here. Otherwise compute the advance
			if ( not spanRecords_ )
				irecord += nBndsVASP_;
			else
				if ( wavefuncDouble_ )
					irecord += nBndsVASP_*(1 + this->num_records_spanned(
							npwSpinKpt_[ik*nspin_+ispin]*sizeof( VASPDprec )*2 ) );
				else
					irecord += nBndsVASP_*(1 + this->num_records_spanned(
							npwSpinKpt_[ik*nspin_+ispin]*sizeof( VASPSprec )*2 ));

		}
}

void ReadVASPWaveFunction::set_up_fourier_max()
{
	//Some helper function to deal with these 3x3 matrices and 3 vectors
	typedef std::vector<double> V;
	auto d_prod = [] (V const & v1, V const & v2 ) {
		return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
	};
	auto sliceBV = [] (V const& matrix3x3, int bv ) {
		return V({matrix3x3[3*0+bv],matrix3x3[3*1+bv],matrix3x3[3*2+bv]});
	};

	//Get the lattice basis and scale to their natural units
	auto A = lattice_.get_latticeMatrix();
	for ( auto &aij : A )
		aij *= lattice_.get_alat();
	auto B = lattice_.get_reciprocal_latticeMatrix();
	for ( auto &bij : B )
		bij *= 2*M_PI/lattice_.get_alat();

	//Check the possible G vectors for plane wave coefficient
	//We compare 4 values - the simple guess below and a modifier
	//where we project on the reciprocal lattice direction.
	double gc = std::sqrt(ecutoff_*energyConverionFactorVASP_);
	std::vector<int> nGmxyz(4, std::floor(gc) + 2 );
	for ( int i = 0 ; i < 3 ; ++i)
		nGmxyz[i] = std::floor(2.0*gc/std::abs(d_prod(sliceBV(A,i),sliceBV(B,i)))
				*sqrt( d_prod(sliceBV(A,i),sliceBV(A,i)))+0.5);
	int nGm = *std::max_element( nGmxyz.begin(),nGmxyz.end() );

	//nGm should now be a lattice vector that is larger than the maximum plane
	//wave in any direction that is below the cutoff in energy.
	//Now we actually walk through them determine the maximum in each direction.
	fourierMax_ = std::vector<int>(3,0);
	V G(3);
	for ( int igx = -nGm; igx <= nGm; ++ igx)
		for ( int igy = -nGm; igy <= nGm; ++ igy)
			for ( int igz = -nGm; igz <= nGm; ++ igz)
			{
				for ( int i = 0 ; i < 3; ++i)
					G[i] = B[i*3+0]*igx+B[i*3+1]*igy+B[i*3+2]*igz;
				if ( d_prod(G,G)/energyConverionFactorVASP_ < ecutoff_ )
				{
					fourierMax_[0] = std::max(fourierMax_[0],std::abs(igx));
					fourierMax_[1] = std::max(fourierMax_[1],std::abs(igy));
					fourierMax_[2] = std::max(fourierMax_[2],std::abs(igz));
				}
			}

	//We need to hold + and - direction (and zero). Also G+k can exceed this
	// number in each direction by up to one so we compute
	for ( auto & nGxi: fourierMax_ )
		nGxi = 2*(nGxi+1)+2;
}

void ReadVASPWaveFunction::read_wavefunction(
		std::vector<int> const& kptindices,
		std::vector<int> const & bandIndices,
		std::vector< std::complex<float> > & wfctData,
		std::vector< int > & npwPerKpt) const
{
	if ( not wavecarfile_.is_open() )
		throw std::logic_error( "Can only read wavefunction from an open file " );

	std::vector<double> B = lattice_.get_reciprocal_latticeMatrix();
	for ( auto &bij : B )
		bij *= 2*M_PI/lattice_.get_alat();

	int numKIrred = static_cast<int>(kpoints_.size())/3;
	int Nk = static_cast<int>(kptindices.size());
	int Nb = static_cast<int>(bandIndices.size());
	assert( Nk <= numKIrred );
	assert( Nb <= nBndsVASP_ );

	//Determine the total amount of storage and the locator for the
	//beginning of each set of plane waves
	int pwPerBandAndSpin = 0;
	if ( int(npwPerKpt.size()) != Nk )
		npwPerKpt.resize(Nk);
	std::vector<int> bandAndKptToPos(nspin_*Nk,0);
	for ( int i = 0 ; i < Nk; ++i)
		for ( int is = 0 ; is < nspin_; ++is)
		{
			int ik = kptindices[i];
			assert( (ik < numKIrred) && ( ik >= 0 ) );
			npwPerKpt[i] = npwSpinKpt_[ik*nspin_+is];
			bandAndKptToPos[i*nspin_+is] = pwPerBandAndSpin*Nb;
			pwPerBandAndSpin += npwSpinKpt_[ik*nspin_+is];
		}
	wfctData.resize( pwPerBandAndSpin*Nb*nspin_ );

	int maxSizePW = *std::max_element(npwPerKpt.begin(),npwPerKpt.end());
	std::vector<char> buffer( maxSizePW*sizeof(VASPDprec)*2 );
	for ( int is = 0 ; is < nspin_; ++is)
		for ( int i = 0 ; i < Nk; ++i)
		{
			int ik = kptindices[i];
			assert( (ik < numKIrred) && ( ik >= 0 ) );

			//locate the base record for this k point and spin
			std::size_t wfctrecordStart = kptSpinPosToFile_[ik*nspin_+is]
							   + this->num_records_spanned( (4+nBndsVASP_*3)*sizeof(VASPDprec) );

			//calculate how many records we span per wavefunction records
			std::size_t byteThisWavefunction;
			if ( wavefuncDouble_ )
				byteThisWavefunction = npwSpinKpt_[ik*nspin_+is]*sizeof(VASPDprec)*2;
			else
				byteThisWavefunction = npwSpinKpt_[ik*nspin_+is]*sizeof(VASPSprec)*2;
			int recPerWfct = 1 + this->num_records_spanned( byteThisWavefunction );;

			for (int j = 0 ; j < Nb; ++j)
			{
				assert( (bandIndices[j] < nBndsVASP_) && ( bandIndices[j] >= 0 ) );

				//here we compute the record locator offset from the ik and spin start
				int rec = wfctrecordStart + recPerWfct*bandIndices[j];
				assert( rec * recl_ + byteThisWavefunction <= total_size_) ;

				wavecarfile_.seekg( rec * recl_ );
				wavecarfile_.read( &buffer[0], byteThisWavefunction  );

				//Now interpret the buffer either as double of float and copy it to the array with the
				//wavefunction data.
				typedef std::complex<VASPDprec> const * VPD;
				typedef std::complex<VASPSprec> const * VPS;

				for ( int ipw = 0 ; ipw < npwSpinKpt_[ik*nspin_+is]; ++ipw)
					wfctData[bandAndKptToPos[i*nspin_+is]+j*npwSpinKpt_[ik*nspin_+is]+ipw ] = wavefuncDouble_ ?
								static_cast< std::complex<float> >(reinterpret_cast<VPD>( &buffer[0] )[ipw] ) :
								static_cast< std::complex<float> >(reinterpret_cast<VPS>( &buffer[0] )[ipw] ) ;
			}
		}
}

void
ReadVASPWaveFunction::compute_fourier_map(
		std::vector<double> const& kptCoords,
		std::vector< std::vector<int> > & fftMapPerK) const
{
	assert( kptCoords.size() % 3 == 0 );
	int Nk = static_cast<int>(kptCoords.size()/3);
	if ( static_cast<int>(fftMapPerK.size()) != Nk )
		fftMapPerK=std::vector< std::vector<int> >(Nk);

	//here we generate the fourier map
	std::vector<double> kPlusG = { 0.0, 0.0, 0.0 };
	std::vector<double> B = lattice_.get_reciprocal_latticeMatrix();
	for ( auto &bij : B )
		bij *= 2*M_PI/lattice_.get_alat();
	for ( int is = 0 ; is < nspin_; ++is)
		for ( int ik = 0 ; ik < Nk; ++ik)
		{
			fftMapPerK[ik] = std::vector<int>(fourierMax_[0]*fourierMax_[1]*fourierMax_[2]*3);
			int ng = 0;
			for ( int igz = 0 ; igz < fourierMax_[2]; ++igz)
			{
				int igzf = igz < fourierMax_[2]/2 ? igz : igz - fourierMax_[2];
				for ( int igy = 0 ; igy < fourierMax_[1]; ++igy)
				{
					int igyf = igy < fourierMax_[1]/2 ? igy : igy - fourierMax_[1];
					for ( int igx = 0 ; igx < fourierMax_[0]; ++igx)
					{
						int igxf = igx < fourierMax_[0]/2 ? igx : igx - fourierMax_[0];
						for ( int xi = 0 ; xi < 3; ++xi)
							kPlusG[xi] = B[xi*3+0]*(kptCoords[ik*3+0]+igxf)
										+B[xi*3+1]*(kptCoords[ik*3+1]+igyf)
										+B[xi*3+2]*(kptCoords[ik*3+2]+igzf);
						if ( (kPlusG[0]*kPlusG[0] + kPlusG[1]*kPlusG[1]
							  + kPlusG[2]*kPlusG[2])/energyConverionFactorVASP_ < ecutoff_ )
						{
							fftMapPerK[ik][ng*3+0] = igx;
							fftMapPerK[ik][ng*3+1] = igy;
							fftMapPerK[ik][ng*3+2] = igz;
							ng++;
						}
					}
				}
			}
			fftMapPerK[ik].resize(ng*3);
		}
}

std::vector<int>
ReadVASPWaveFunction::get_fft_max_dims() const
{
	return fourierMax_;
}

int ReadVASPWaveFunction::num_records_spanned(std::size_t bytesToRead) const
{
	return (static_cast<std::int64_t>(bytesToRead)-1)/static_cast<std::int64_t>(recl_);
}

std::vector<double> const& ReadVASPWaveFunction::get_energies() const
{
	return energies_;
}

int ReadVASPWaveFunction::get_num_bands() const
{
	return nBndsVASP_;
}

int ReadVASPWaveFunction::get_num_spins() const
{
	return nspin_;
}

int ReadVASPWaveFunction::get_num_kpts() const
{
	return static_cast<int>(kpoints_.size())/3;
}

std::vector<double> const & ReadVASPWaveFunction::get_k_points() const
{
	return kpoints_;
}

} /* namespace IOMethods */
} /* namespace elephon */
