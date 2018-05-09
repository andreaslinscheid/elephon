/*	This file FrozenCore.h is part of elephon.
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
 *  Created on: Apr 27, 2018
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_ATOMICSITE_FROZENCORE_H_
#define ELEPHON_ATOMICSITE_FROZENCORE_H_

#include "AtomicSite/RadialGrid.h"
#include <vector>

namespace elephon
{
namespace AtomicSite
{

/**
 * Class to represent the frozen core contribution including the point charge of the protons and the core electronic charge.
 */
class FrozenCore
{
public:

	/**
	 *	Initialize the FrozenCore class.
	 *
	 * @param[in] corePointCharge				Number of protons in the core of the atom
	 * @param[in] electronicFrozenCoreCharge	The charge density of the electronic frozen core system on the \p radialGrid times r^2 with
	 * 											r in units of Angstrom.
	 * @param[in] radialGrid					The radial grid.
	 */
	void initialize(double corePointCharge,
			std::vector<double> electronicFrozenCoreCharge,
			RadialGrid radialGrid);
	/**
	 * Add the core point charge to the spherical potential.
	 *
	 * The potential will be given in terms of an expansion in spherical harmonics, but is spherical, so it will be divided by \f$\sqrt{4\pi}\f$
	 *
	 * @param[in] sphericalPotentialToBeAddedToBegin	Add the -Z/r potential to the container pointed to by sphericalPotentialToBeAddedToBegin.
	 * 													The target should be expansion coefficients of the spherical harmonic l=0,m=0 and be measured in units
	 * 													of eV.
	 * @param[in] sphericalPotentialToBeAddedToEnd		End of the range of the potential to be added to. The range must be equal to the internal radialGrid.
	 */
	template<class iterator>
	void add_core_potential(
			iterator sphericalPotentialToBeAddedToBegin,
			iterator sphericalPotentialToBeAddedToEnd) const;

	/**
	 * Add the Hartree potential of the frozen core charge to the spherical potential.
	 *
	 * The potential will be given in terms of an expansion in spherical harmonics, but is spherical, so it will be divided by \f$\sqrt{4\pi}\f$
	 *
	 * @param[in] sphericalPotentialToBeAddedToBegin	Add the Hartree potential to the container pointed to by sphericalPotentialToBeAddedToBegin.
	 * 													The target should be expansion coefficients of the spherical harmonic l=0,m=0 and be measured in units
	 * 													of eV.
	 * @param[in] sphericalPotentialToBeAddedToEnd		End of the range of the potential to be added to. The range must be equal to the internal radialGrid.
	 */
	template<class iterator>
	void add_core_hartree_potential(
			iterator sphericalPotentialToBeAddedToBegin,
			iterator sphericalPotentialToBeAddedToEnd) const;

	/**
	 * Add the derivative of Hartree potential of the frozen core charge and the core point charge to the displacement potential.
	 *
	 * The displacement potential will be given in terms of an expansion in spherical harmonics. Since the bare potential is spherically symmetryic
	 * the displacement will in the l=1 channels. The projection  so it will be divided by \f$\sqrt{4\pi}\f$
	 *
	 * @param[in] direction						Specify the displacement direction in Cartesian coordinates w.r.p.t. the coordiante system of the radial grid.
	 * 											It may or may not be normalized to 1.0. If e.g. all 3 channels x,y, and z are needed, simply pass direction (1,1,1).
	 * @param[in] potentialToBeAddedToBegin		Add the frozen core potential to the container pointed to by potentialToBeAddedToBegin.
	 * 											The target should be expansion coefficients of the spherical harmonic in the
	 * 											layout radial grid as the fasted running index, then the layout specified by
	 * 											Auxilliary::memlayout::angular_momentum_layout().
	 * 											The result will be added to the l=1 channels. Input must therefore include those channels.
	 * @param[in] potentialToBeAddedToEnd		End of the range of the potential to be added to.
	 */
	template<class iterator>
	void add_core_hartree_analytic_displacement(
			std::array<double,3> direction,
			iterator potentialToBeAddedToBegin,
			iterator potentialToBeAddedToEnd) const;

	/**
	 * Compute the total charge contained in this object.
	 * @return The charge in units of the electorn elementary charge.
	 */
	double total_electronic_charge() const;

	/**
	 * Inspect the radial grid of this object.
	 * @return a constant reference to the internal radial grid.
	 */
	RadialGrid const & get_radial_grid() const;
private:

	double corePointCharge_ = 0.0;

	std::vector<double> electronicFrozenCoreCharge_;

	RadialGrid rgrid_;

	mutable std::vector<double> coreHartreeZDisplacementDataBuffer_;

	void compute_core_hartree_analytic_z_displacement() const;


	void set_hartree_displ(int ir, double radialIntegralCharge) const;

	void add_contribution_prev_two_intervals(int ir, double & integral, std::vector<double> const & y) const;
};

} /* namespace AtomicSite */
} /* namespace elephon */

#include "AtomicSite/FrozenCore.hpp"
#endif /* ELEPHON_ATOMICSITE_FROZENCORE_H_ */
