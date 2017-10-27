/*	This file AlphaSquaredF.h is part of elephon.
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
 *  Created on: Sep 26, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_PHONONSTRUCTURE_ALPHASQUAREDF_H_
#define ELEPHON_PHONONSTRUCTURE_ALPHASQUAREDF_H_

#include "IOMethods/ResourceHandler.h"
#include <memory>

namespace elephon
{
namespace PhononStructure
{

class AlphaSquaredF
{
public:

	/**
	 * This is the main driver for the a2F calculation.
	 *
	 * Computes and sets the a2F coupling function internally.
	 *
	 * @param resourceHandler	The interface on how to obtain all the data need for this task.
	 */
	void compute_a2F( std::shared_ptr<IOMethods::ResourceHandler> resourceHandler );

	void write_a2F_file(std::string const & filename) const;
private:

	std::vector<double> a2F_;

	double freqMin_ = 0;

	double freqMax_ = 0;

	int freqNPts_ = 0;

	void map_freq_grid_slot(
			std::vector<float>::const_iterator begData,
			std::vector<float>::const_iterator begFreq,
			std::vector<float>::const_iterator endFreq);
};

} /* namespace PhononStructure */
} /* namespace elephon */

#endif /* ELEPHON_PHONONSTRUCTURE_ALPHASQUAREDF_H_ */
