/*	This file NestingFunction.h is part of elephon.
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
 *  Created on: Oct 22, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_ELECTRONICSTRUCTURE_NESTINGFUNCTION_H_
#define ELEPHON_ELECTRONICSTRUCTURE_NESTINGFUNCTION_H_

#include "ElectronicStructure/TetrahedraIsosurface.h"
#include <vector>
#include <memory>

namespace elephon
{
namespace ElectronicStructure
{

class NestingFunction
{
public:

	void compute_nesting_function(
			std::shared_ptr<const TetrahedraIsosurface> isoSurface);

	void generate_vtk_volume_data();

	std::shared_ptr<LatticeStructure::DataRegularGrid<float> > get_nesting(int isoEnergyIndex) const;

private:

	std::vector< std::shared_ptr<LatticeStructure::DataRegularGrid<float> > > nesting_;

	std::shared_ptr<const TetrahedraIsosurface> isoSurface_;
};

} /* namespace ElectronicStructure */
} /* namespace elephon */

#endif /* ELEPHON_ELECTRONICSTRUCTURE_NESTINGFUNCTION_H_ */
