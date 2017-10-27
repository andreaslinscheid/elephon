/*	This file BandStructureAnalysis.h is part of elephon.
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
 *  Created on: Sep 7, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_ELECTRONICSTRUCTURE_BANDSTRUCTUREANALYSIS_H_
#define ELEPHON_ELECTRONICSTRUCTURE_BANDSTRUCTUREANALYSIS_H_

#include "ElectronicStructure/ElectronicBands.h"
#include "IOMethods/ResourceHandler.h"
#include <string>
#include <vector>
#include <memory>

namespace elephon
{
namespace ElectronicStructure
{
namespace BandStructureAnalysis
{

/**
 * A type to hold a single band extrema containing a band index
 * and the k point location in the Brillouin zone.
 */
typedef struct {
	int ibnd;
	std::vector<double> k;
	double energy;
} b_extrema;

void write_mass_tensor_file(
		std::string const & filename,
		ElectronicBands const & bands,
		int startBand,
		std::vector<double> const & energyWindow,
		int method);

void do_band_structure_analysis(std::shared_ptr<IOMethods::ResourceHandler> loader);

} /* namespace BandStructureAnalysis */
} /* namespace ElectronicStructure */
} /* namespace elephon */

#endif /* ELEPHON_ELECTRONICSTRUCTURE_BANDSTRUCTUREANALYSIS_H_ */
