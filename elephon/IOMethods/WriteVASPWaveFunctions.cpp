/*	This file WriteVASPWaveFunctions.cpp is part of elephon.
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
 *  Created on: Aug 10, 2017
 *      Author: A. Linscheid
 */

#include "IOMethods/WriteVASPWaveFunctions.h"
#include "IOMethods/ReadVASPWaveFunction.h"

namespace elephon
{
namespace IOMethods
{

void write_VASP_wavefunctions(
		std::string const & filename,
		double eCutPW,
		LatticeStructure::LatticeModule const & lattice,
		ElectronicStructure::ElectronicBands const & bands,
		std::vector<std::vector<std::complex<float>>> const & wavefunctions,
		std::vector<std::vector<int>> const & fftMap )
{	// We create a file that has the format of the VASP wavecar.
	// VASP writes double precision (Fortran selected_real_kind(10))
	// or single precision (Fortran selected_real_kind(5)) dependent on the version of the code.
	typedef double VASPDP;
	typedef float VASPSP;

	const int nkpIrred = bands.get_grid().get_np_irred();
	LatticeStructure::RegularSymmetricGrid const & kgrid = bands.get_grid();
	int nBnd = bands.get_nBnd();

	std::ofstream file( filename.c_str() , std::ios_base::binary );
	if ( not file.good() )
		throw std::runtime_error( std::string("Error opening file ")+filename);

	//Header done. Now we write the wavefunctions, first we compute the VASP fft maps
	std::vector<std::vector<int>> fftMapVASP(nkpIrred);
	std::vector<int> npwK(nkpIrred);

	std::vector<int> fourierMax;
	ReadVASPWaveFunction reader;
	reader.compute_fourier_max(
			eCutPW,
			lattice,
			fourierMax);

	std::vector<double> kvectors(3*nkpIrred);
	for ( int ikir = 0 ; ikir < nkpIrred; ++ikir)
	{
		//Create the plane wave set at this k point
		int ikred = kgrid.get_maps_irreducible_to_reducible()[ikir][ kgrid.get_symmetry().get_identity_index() ];
		auto k = kgrid.get_vector_direct(ikred);
		// use the vasp convention for k points
		for ( int i = 0 ; i < 3; ++i)
			kvectors[ikir*3+i] =  k[i] < -0.5 + 0.5*kgrid.get_grid_prec() ?
					k[i] + 1.0 - 0.5*kgrid.get_grid_prec() : k[i];
	}
	reader.compute_fourier_map(
			kvectors,
			fftMapVASP,
			kgrid.get_grid_prec(),
			bands.get_nspin(),
			fourierMax,
			eCutPW,
			lattice);

	//TODO this does not work for span records in VASP
	int reclength = std::max(3+9, nBnd*3+4)*sizeof(VASPDP);//In bytes
	for (int ik = 0 ; ik < nkpIrred; ++ik)
	{
		int wfctBlock = 2*(fftMapVASP[ik].size()/3)*sizeof(VASPSP);
		reclength = reclength < wfctBlock ? wfctBlock : reclength;
	}

	//Compute the total number of records in the file as
	//the header plus per k point one record for bands ect and for each
	//band one record for the wavefunction data
	int numRecs = 2 + nkpIrred*(1+bands.get_nBnd());
	std::vector<char> binaryBuffer( numRecs*reclength );
	VASPDP * buffAsFloat = reinterpret_cast<VASPDP * >( &binaryBuffer[0] );

	//Write header
	//First record
	buffAsFloat[0] = reclength;
	buffAsFloat[1] = 1; // #spin
	buffAsFloat[2] = 45200;//Version tag

	//second record
	buffAsFloat = reinterpret_cast<VASPDP * >( &binaryBuffer[reclength] );
	buffAsFloat[0] = nkpIrred;
	buffAsFloat[1] = bands.get_nBnd();
	buffAsFloat[2] = eCutPW;
	for ( int i = 0 ; i < 3; ++i )
		for ( int j = 0 ; j < 3; ++j )
			buffAsFloat[3+i*3+j] = lattice.get_alat()*lattice.get_latticeMatrix()[i*3+j];

	int irecord = 2;
	for ( int ikir = 0 ; ikir < nkpIrred; ++ikir)
	{
		int npw = fftMapVASP[ikir].size()/3;

		if ( npw != wavefunctions[ikir].size() )
			throw std::runtime_error("Number of plane waves passed does not match the VASP number");
		if ( fftMapVASP[ikir].size() != fftMap[ikir].size() )
			throw std::runtime_error("Number of plane waves passed does not match");

		//Create the plane wave set at this k point
		//Helper G Vector struct that follows the VASP ordering
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
		std::map< GVector, int > pwconversion;
		auto hint = pwconversion.end();
		for ( int ipw = 0 ; ipw < npw; ++ipw )
			hint = pwconversion.insert(hint,
					std::make_pair(GVector(fftMapVASP[ikir][ipw*3 + 0], fftMapVASP[ikir][ipw*3 + 1], fftMapVASP[ikir][ipw*3 + 2]),
							ipw));

		buffAsFloat = reinterpret_cast<VASPDP * >( &binaryBuffer[irecord*reclength] );
		buffAsFloat[0] = npw;
		buffAsFloat[1] = kvectors[ikir*3 + 0];
		buffAsFloat[2] = kvectors[ikir*3 + 1];
		buffAsFloat[3] = kvectors[ikir*3 + 2];
		for (int ibnd = 0 ; ibnd < nBnd; ++ibnd)
		{
			buffAsFloat[4+ibnd*3+0] = bands(ikir, ibnd);
			buffAsFloat[4+ibnd*3+1] = 0.0;// complex part of the eigenvalue ... ?!
			buffAsFloat[4+ibnd*3+2] = 1.0;// occupation is neglected
		}

		for (int ibnd = 0 ; ibnd < nBnd; ++ibnd)
		{
			VASPSP * buffAsSingle = reinterpret_cast<VASPSP * >( &binaryBuffer[(irecord + ibnd + 1)*reclength] );
			for ( int ipw = 0 ; ipw < npw; ++ipw)
			{
				auto it = pwconversion.find(
						GVector(fftMap[ikir][ipw*3 + 0], fftMap[ikir][ipw*3 + 1], fftMap[ikir][ipw*3 + 2]));
				if ( it == pwconversion.end() )
					throw std::runtime_error("Problem locating the correct plane wave coefficient");
				buffAsSingle[(it->second)*2+0] = std::real(wavefunctions[ikir][ibnd*npw+ipw]);
				buffAsSingle[(it->second)*2+1] = std::imag(wavefunctions[ikir][ibnd*npw+ipw]);
			}
		}
		irecord += 1 + nBnd;
	}
	file.write(binaryBuffer.data(),binaryBuffer.size());
	file.close();
}

} /* namespace IOMethods */
} /* namespace elephon */
