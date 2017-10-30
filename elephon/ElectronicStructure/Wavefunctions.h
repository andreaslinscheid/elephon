/*	This file Wavefunctions.h is part of elephon.
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

#ifndef ELEPHON_ELECTRONICSTRUCTURE_WAVEFUNCTIONS_H_
#define ELEPHON_ELECTRONICSTRUCTURE_WAVEFUNCTIONS_H_

#include "IOMethods/ElectronicStructureCodeInterface.h"
#include "LatticeStructure/RegularSymmetricGrid.h"
#include <vector>
#include <complex>
#include <memory>

namespace elephon
{
namespace ElectronicStructure
{

class Wavefunctions
{
public:
	void initialize(
			LatticeStructure::RegularSymmetricGrid kgrid,
			std::shared_ptr<IOMethods::ElectronicStructureCodeInterface> wfctInterface);

	void initialize(
			std::string rootDir,
			std::shared_ptr<IOMethods::ElectronicStructureCodeInterface> wfctInterface);

	void initialize(
			std::string rootDir,
			LatticeStructure::RegularSymmetricGrid kgrid,
			int numBands,
			std::shared_ptr<IOMethods::ElectronicStructureCodeInterface> wfctInterface);

	void generate_reducible_grid_wfcts(
			std::vector<int> const & bndIndices,
			std::vector<int> const & redKptIndices,
			std::vector< std::vector<std::complex<float>> > & wfcts,
			std::vector<int> & npwPerKAndSpin) const;

	int get_num_bands() const;

	LatticeStructure::RegularSymmetricGrid const & get_k_grid() const;

	void compute_Fourier_maps(std::vector<double> const & kvectors,
			std::vector<std::vector<int>> & fftMapsPerK) const;

	void generate_wfcts_at_arbitray_kp(
			std::vector<double> kList,
			std::vector<int> const & bandList,
			std::vector< std::vector< std::complex<float> > > & wfctsArbitrayKp,
			std::vector<std::vector<int>> & fftMapsArbitrayKp) const;

	/**
	 * Obtain the Fourier grid sufficient to represent any wavefunction in the system.
	 *
	 * @return vector with max x, max y and max z reciprocal lattice vectors
	 */
	std::vector<int> get_max_fft_dims() const;
private:

	typedef class KPGVect
	{
	public:
		KPGVect(double prec, double kpG[3]) : gridPrec_(prec), kpG_{kpG[0],kpG[1],kpG[2]} { };
		bool operator< (KPGVect const & other) const
		{
			assert( gridPrec_ == other.gridPrec_ );
			for ( int i = 0 ; i < 3 ; ++i )
				if ( std::abs(kpG_[i]-other.kpG_[i]) > gridPrec_ )
					return kpG_[i] < other.kpG_[i];
			return false;
		}
	private:
		double gridPrec_ = 0;
		double kpG_[3] = {0,0,0};
	} KPlusGVector;

	int nBnd_ = -1;

	std::string rootDir_;

	std::shared_ptr<IOMethods::ElectronicStructureCodeInterface> wfctInterface_;

	LatticeStructure::RegularSymmetricGrid grid_;

	mutable std::vector< std::vector<int> > gSymBuffer_;

	mutable std::vector< std::vector< std::complex<float> > > phaseBuffer_;

	void fill_G_symmetry_buffer(int isym) const;
};

} /* namespace ElectronicStructure */
} /* namespace elephon */

#endif /* ELEPHON_ELECTRONICSTRUCTURE_WAVEFUNCTIONS_H_ */
