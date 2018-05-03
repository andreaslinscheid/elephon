/*	This file AtomSiteData.h is part of elephon.
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
 *  Created on: Feb 2, 2018
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_ATOMICSITE_ATOMSITEDATA_H_
#define ELEPHON_ATOMICSITE_ATOMSITEDATA_H_

#include "LatticeStructure/Atom.h"
#include "AtomicSite/SphericalHarmonicExpansion.h"
#include "AtomicSite/FrozenCore.h"

namespace elephon
{
namespace AtomicSite
{

/**
 * A collection of an atomic object and the spherical harmonic expansion of data around it.
 */
class AtomSiteData
{
public:

	/**
	 * Initialize this container for atomic site data.
	 *
	 * @param[in] a
	 * @param[in] dataPotential
	 * @param[in] coreData
	 */
	void initialize(
			LatticeStructure::Atom a,
			SphericalHarmonicExpansion dataPotential,
			FrozenCore coreData);

	LatticeStructure::Atom const & get_atom() const;

	SphericalHarmonicExpansion const & get_potential_data() const;

	SphericalHarmonicExpansion & edit_potential_data();

	FrozenCore const & get_frozen_core_data() const;
private:

	LatticeStructure::Atom a_;

	SphericalHarmonicExpansion dataPotential_;

	FrozenCore coreData_;
};

} /* namespace AtomicSite */
} /* namespace elephon */

#endif /* ELEPHON_ATOMICSITE_ATOMSITEDATA_H_ */
