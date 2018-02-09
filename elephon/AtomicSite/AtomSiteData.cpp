/*	This file AtomSiteData.cpp is part of elephon.
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

#include "AtomicSite/AtomSiteData.h"

namespace elephon
{
namespace AtomicSite
{

void
AtomSiteData::initialize(
		LatticeStructure::Atom a,
		SphericalHarmonicExpansion data)
{
	a_ = std::move(a);
	data_ = std::move(data);
}

LatticeStructure::Atom const &
AtomSiteData::get_atom() const
{
	return a_;
}

SphericalHarmonicExpansion const &
AtomSiteData::get_data() const
{
	return data_;
}

SphericalHarmonicExpansion &
AtomSiteData::edit_data()
{
	return data_;
}

} /* namespace AtomicSite */
} /* namespace elephon */
