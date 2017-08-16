/*	This file ReadVASPWaveFunction.h is part of elephon.
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

#ifndef ELEPHON_IOMETHODS_READVASPWAVEFUNCTION_H_
#define ELEPHON_IOMETHODS_READVASPWAVEFUNCTION_H_

#include "LatticeStructure/LatticeModule.h"
#include <cstdint>
#include <string>
#include <complex>
#include <vector>
#include <fstream>

namespace elephon
{
namespace IOMethods
{

class ReadVASPWaveFunction
{

public:

	ReadVASPWaveFunction();

	~ReadVASPWaveFunction();

	void prepare_wavecar(std::string filename);

	void read_wavefunction(
			std::vector<int> const& kptindices,
			std::vector<int> const & bandIndices,
			std::vector< std::vector< std::complex<float> > > & wfctData,
			std::vector<int> & npwPerKpt ) const;

	void compute_fourier_map(
			std::vector<double> kptCoords,
			std::vector< std::vector<int> > & fftMapPerK,
			double vaspGridPrec) const;

	void compute_fourier_map(
			std::vector<double> kptCoords,
			std::vector< std::vector<int> > & fftMapPerK,
			double vaspGridPrec,
			int nspin,
			std::vector<int> const & fourierMax,
			double ecutoff,
			LatticeStructure::LatticeModule const & lattice) const;

	std::vector<int> const & get_fft_max_dims() const;

	void compute_fourier_max(
			double ecutoff,
			LatticeStructure::LatticeModule const & lattice,
			std::vector<int> & fourierMax) const;

	std::vector<double> const & get_energies() const;

	int get_num_bands() const;

	int get_num_spins() const;

	int get_num_kpts() const;

	std::vector<double> const & get_k_points() const;

	std::string const & get_filename() const;
private:

	typedef double VASPDprec;

	typedef float VASPSprec;

	std::size_t recl_ = 0;

	std::size_t total_size_ = 0;

	int nBndsVASP_ = 0;

	double ecutoff_ = 0;

	int nspin_ = 1;

	bool wavefuncDouble_ = false;

	bool spanRecords_ = false;

	/// vector of k-point coordinates as (x1,y1,z1,x2,y2,...)
	///
	/// NOTE: elephon stores k points in the 'left major' 1.BZ so that k=(0.5 0.0 -0.5) will
	/// be mapped to k_1BZ=(-0.5 0.0 -0.5) while VASP apparently uses k_1BZ=(0.5 0.0 0.5).
	/// We have to keep track of this while reading in data and output only data in the
	/// elephon convention! This affects also the Fourier mapping since since elephon's G vectors
	/// for k points on the border are off by possibly 1 in each direction.
	std::vector<double> kpoints_;

	std::vector<int> fourierMax_;

	elephon::LatticeStructure::LatticeModule lattice_;

	//The file is kept open and closed in the destructor
	mutable std::fstream wavecarfile_;

	/**
	 * For each k point index in the regular
	 */
	std::vector< std::vector<std::int64_t> > byteLocationMap_;

	std::vector<int> npwSpinKpt_;

	std::string filename_;

	std::vector<double> energies_;

	//Element (ik,ispin) = ik*nspin_+ispin returns the record number where to find the beginning of the
	// data for k point 'ik' and spin 'ispin'. The wavefunctions start one (or more) later.
	std::vector<std::size_t> kptSpinPosToFile_;

	const double energyConverionFactorVASP_ = 0.262465831;

	void set_up_fourier_max();

	int num_records_spanned(std::size_t bytesToRead) const;
};

} /* namespace IOMethods */
} /* namespace elephon */

#endif /* ELEPHON_IOMETHODS_READVASPWAVEFUNCTION_H_ */
