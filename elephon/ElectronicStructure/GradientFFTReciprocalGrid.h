/*	This file GradientFFTReciprocalGrid.h is part of elephon.
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
 *  Created on: Apr 29, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_ELECTRONICSTRUCTURE_GRADIENTFFTRECIPROCALGRID_H_
#define ELEPHON_ELECTRONICSTRUCTURE_GRADIENTFFTRECIPROCALGRID_H_

#include "LatticeStructure/LatticeModule.h"
#include <vector>
#include <cstdlib>

namespace elephon
{
namespace ElectronicStructure
{

class GradientFFTReciprocalGrid {
public:
	GradientFFTReciprocalGrid();

	void compute_gradient(
			std::vector<int> grid,
			LatticeStructure::LatticeModule const& lattice,
			int nDataBlock,
			std::vector<double> const& dataOnGrid);

	void copy_data(std::vector<int> const& conseqGridIndices,
			std::vector<int> bands,
			std::vector<double> & dataAtIndices) const;

	std::vector<double> const& get_data() const;
private:
	std::vector<int> grid_;

	int nBlockData_ = 0;

	std::vector<double> gradientDataOnGrid_;

};

} /* namespace ElectronicStructure */
} /* namespace elephon */

#endif /* ELEPHON_ELECTRONICSTRUCTURE_GRADIENTFFTRECIPROCALGRID_H_ */
