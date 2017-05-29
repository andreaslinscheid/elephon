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

#include "LatticeStructure/RegularGrid.h"
#include "IOMethods/ElectronicStructureCodeInterface.h"
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
			double gridPrec,
			std::string rootDir,
			std::shared_ptr<IOMethods::ElectronicStructureCodeInterface> wfctInterface);

	void generate_reducible_grid_wfcts(
			std::vector<int> const & bndIndices,
			std::vector<int> const & redKptIndices,
			std::vector< std::complex<float> > & wfcts,
			std::vector<int> & npwPerKAndSpin) const;
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

	int nBnd_;

	std::string rootDir_;

	std::vector< std::vector<int> > fourierMap_;

	std::shared_ptr<IOMethods::ElectronicStructureCodeInterface> wfctInterface_;

	LatticeStructure::RegularGrid grid_;

	mutable std::vector< std::vector<int> > gSymBuffer_;

	mutable std::vector< std::vector< std::complex<float> > > phaseBuffer_;

	void fill_G_symmetry_buffer(int isym) const;
};

} /* namespace ElectronicStructure */
} /* namespace elephon */

#endif /* ELEPHON_ELECTRONICSTRUCTURE_WAVEFUNCTIONS_H_ */
