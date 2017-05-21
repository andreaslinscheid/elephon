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
	void initialize(
			std::vector<double> const & kpoints,
			int numBands,
			std::shared_ptr<IOMethods::ElectronicStructureCodeInterface> wfctInterface,
			LatticeStructure::RegularGrid grid);

	void generate_reducible_grid_wfcts(
			std::vector<int> const & bndIndices,
			std::vector<int> const & redKptIndices,
			std::vector< std::complex<float> > & wfcts) const;

	void get_Fourier_maps(
			std::vector<int> const & redKptIndices,
			std::vector< std::vector<int> > & fftMap) const;
private:

	int nBnd_;

	std::vector< std::vector<int> > fourierMap_;

	std::shared_ptr<IOMethods::ElectronicStructureCodeInterface> wfctInterface_;
};

} /* namespace ElectronicStructure */
} /* namespace elephon */

#endif /* ELEPHON_ELECTRONICSTRUCTURE_WAVEFUNCTIONS_H_ */
