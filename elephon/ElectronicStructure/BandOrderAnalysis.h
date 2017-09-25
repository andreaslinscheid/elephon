/*	This file BandOrderAnalysis.h is part of elephon.
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
 *  Created on: Sep 11, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_ELECTRONICSTRUCTURE_BANDORDERANALYSIS_H_
#define ELEPHON_ELECTRONICSTRUCTURE_BANDORDERANALYSIS_H_

#include "ElectronicStructure/Wavefunctions.h"
#include <vector>

namespace elephon
{
namespace ElectronicStructure
{

class BandOrderAnalysis
{
public:

	void compute_band_order_overlap(
			ElectronicBands const & bands,
			Wavefunctions const & wfct);

	int operator() (int ikir, int ib) const;
private:

	const double energyGradientCutoff_ = 15.0; // eV * \AA

	int nB_ = 0;

	std::vector<int> bandOrder_;

	void compute_irred_BZ_paths(
			LatticeStructure::RegularSymmetricGrid const& kgrid,
			std::vector<std::vector<int>> & paths) const;

	void compute_band_overlap_matrix(
			Wavefunctions const & wfct,
			int ikr1, std::vector<int> const & bands1,
			int ikr2, std::vector<int> const & bands2,
			std::vector<float> & overlaps) const;
};

} /* namespace ElectronicStructure */
} /* namespace elephon */

#endif /* ELEPHON_ELECTRONICSTRUCTURE_BANDORDERANALYSIS_H_ */
