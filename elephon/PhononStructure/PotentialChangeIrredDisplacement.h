/*	This file PotentialChangeIrredDisplacement.h is part of elephon.
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
 *  Created on: Jan 17, 2018
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_PHONONSTRUCTURE_POTENTIALCHANGEIRREDDISPLACEMENT_H_
#define ELEPHON_PHONONSTRUCTURE_POTENTIALCHANGEIRREDDISPLACEMENT_H_

#include "Auxillary/AlignedVector.h"
#include "LatticeStructure/DataRegularAndRadialGrid.h"
#include <memory>
#include <vector>

namespace elephon
{
// forward declares
namespace LatticeStructure { template<typename T> class DataRegularAndRadialGrid; };
namespace LatticeStructure { class RegularBareGrid; };
namespace LatticeStructure { class UnitCell; };
namespace LatticeStructure { class AtomDisplacement; };
namespace LatticeStructure { class PrimitiveToSupercellConnection; };

namespace PhononStructure
{

/**
 * Handles the determination of the potential difference induced by the irreducible displacements.
 *
 * From the raw input data, we determine for a irreducible atomic displacement the change in potential
 * on the unperturbed grid.
 * The initializing method takes data for the perturbed and the unperturbed atomic configuration. Since the regular grid
 * data part is on the same grid, the resulting finite difference is simply
 *
 * \f{eqnarray*}{
 * \Delta v_{{\rm {\scriptscriptstyle scf}}}({\bf r}_{i})
 * 		& = &
 * 	v_{{\rm {\scriptscriptstyle scf}}}[\{{\bf \tau}_{\kappa{\bf R}}^{0}\}]({\bf r}_{i})
 * 			-v_{{\rm {\scriptscriptstyle scf}}}[\ldots,{\bf \tau}_{\kappa_{1}{\bf R}_{1}}^{0},
 * 			  {\bf \tau}_{\kappa\alpha{\bf R}}^{0}+\Delta{\bf \tau}_{\kappa\alpha{\bf R}},{\bf \tau}_{\kappa_{2}{\bf R}_{2}}^{0},\ldots]({\bf r}_{i})
 * \f}
 *
 * where \f$ {\bf r}_{i} \f$ is a regular grid vector. Now for the radial data, we are facing the problem that one of sites is shifted. For the radial part
 * we have the following equation
 *
 * \f{eqnarray*}{
 * \Delta v_{{\rm {\scriptscriptstyle scf}}}^{\kappa_{0}}({\bf r})
 * 		& = &
 * 	\sum_{l=0}^{l_{{\rm Max}}}\sum_{m=-l}^{l}\bigl(v_{l,m}(\vert{\bf r}-\tau_{\kappa{\bf r}}^{0}\vert)Y_{l}^{m}(\frac{{\bf r}-\tau_{\kappa{\bf r}}^{0}}{\vert{\bf r}-\tau_{\kappa{\bf r}}^{0}\vert}) \\
 * 		&   &
 * 	 -v_{l,m}(\vert{\bf r}-\tau_{\kappa\alpha{\bf r}}^{0}-\Delta\tau_{\kappa\alpha{\bf r}}\vert)Y_{l}^{m}(\frac{{\bf r}-\tau_{\kappa\alpha{\bf r}}^{0}-\Delta
 * 	 	\tau_{\kappa\alpha{\bf r}}}{\vert{\bf r}-\tau_{\kappa\alpha{\bf r}}^{0}-\Delta\tau_{\kappa\alpha{\bf r}}\vert})\bigr)
 * \f}
 *
 * Here, we take the approach to compute the second term in the above equation, expanded around \f$ {\bf r}-\tau_{\kappa\alpha{\bf r}}^{0}\f$,
 * and fit different expansion coefficients \f$ v_{l,m}^{\prime}(\vert {\bf r}\vert)\f$ at the center of the original atom.
 * This is done using AtomicSite::SphericalHarmonicExpansion::interpolate() following by a AtomicSite::SphericalHarmonicExpansion::fit_to_data().
 * Then, we can compute the finite difference in terms of the expansion coefficients.
 * For the non-displaced sites, we simply subtract the data values on the same grid
 *
 * \f{eqnarray*}{
 * v_{l,m}^{\Delta\kappa}(\vert{\bf r}-\tau_{\kappa{\bf r}}^{0}\vert)
 * 		& = &
 * 	v_{l,m}^{(0)}(\vert{\bf r}-\tau_{\kappa{\bf r}}^{0}\vert)-v_{l,m}^{(1)}(\vert{\bf r}-\tau_{\kappa{\bf r}}^{0}\vert)
 * \f}
 */
class PotentialChangeIrredDisplacement
{
public:

	/**
	 * Compute the difference potential due to a finite displacement.
	 *
	 * @throw a logic error if displ references an atom that is not in the primitive unit cell
	 * 			or radialDisplacedPotential contains two atoms that overlap.
	 *
	 * @param[in] displ
	 * @param[in] regularGridGroundStatePotential
	 * @param[in] radialGroundStatePotential
	 * @param[in] regularGridDisplacedPotential
	 * @param[in] radialDisplacedPotential
	 * @param[in] unitcellGrid
	 * @param[in] supercellGrid
	 */
	void initialize(
			LatticeStructure::AtomDisplacement const & atomDispl,
			LatticeStructure::DataRegularAndRadialGrid<double> const & groundStatePotential,
			LatticeStructure::DataRegularAndRadialGrid<double> const & displacedPotential,
			std::shared_ptr<const LatticeStructure::RegularBareGrid> unitcellGrid,
			std::shared_ptr<const LatticeStructure::RegularBareGrid> supercellGrid,
			std::shared_ptr<const LatticeStructure::PrimitiveToSupercellConnection> scPrimMap );

	/**
	 * Get number of radial data elements stored in the object.
	 * @return	the number of radial data elements
	 */
	int get_max_num_radial_elements() const;

	/**
	 * Get number of angular momentum channels of the data stored in the object.
	 * @return	the number of angular momentum channels
	 */
	int get_max_num_angular_moment_channels() const;

	/**
	 * Re-shuffel the internal data according to the symmetry operation sop.
	 *
	 * @param[in] sop	Representation of the symmetry operation that will be applied.
	 */
	void transform(symmetry::SymmetryOperation const & sop);

	Auxillary::alignedvector::DV::const_iterator begin_regular_data() const;

	Auxillary::alignedvector::DV::const_iterator end_regular_data() const;

	Auxillary::alignedvector::ZV::const_iterator begin_radial_data(int atomIndex) const;

	Auxillary::alignedvector::ZV::const_iterator end_radial_data(int atomIndex) const;

private:

	LatticeStructure::DataRegularAndRadialGrid<double> data_;
};

} /* namespace PhononStructure */
} /* namespace elephon */

#endif /* ELEPHON_PHONONSTRUCTURE_POTENTIALCHANGEIRREDDISPLACEMENT_H_ */
