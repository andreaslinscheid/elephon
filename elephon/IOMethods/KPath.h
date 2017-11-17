/*	This file KPath.h is part of elephon.
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
 *  Created on: Nov 1, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_IOMETHODS_KPATH_H_
#define ELEPHON_IOMETHODS_KPATH_H_

#include "Auxillary/AlignedVector.h"
#include <string>
#include <vector>

namespace elephon
{
namespace IOMethods
{

class KPath
{
public:

	void read_kpath_file(std::string const & filename);

	void produce_gnuplot_script_spectral_function(
			std::string filenameScript,
			std::string filenameData,
			std::string quantityLabel,
			std::vector<double> const & frequencies,
			std::vector<double> const & data) const;

	void produce_gnuplot_script_stable_particle(
			std::string filenameScript,
			std::string filenameData,
			std::string quantityLabel,
			Auxillary::alignedvector::DV const & data,
			int numBnds,
			std::pair<double,double> energyRange) const;

	void save_bands_path_stable_particle(
			std::string const & filename,
			Auxillary::alignedvector::DV const & data,
			int numBnds) const;

	std::vector<double> const & get_k_points() const;
private:

	std::vector<double> kpts_;

	//Note that we don't need labels_ on processors other than the ioproc
	std::vector<std::pair<int,std::string> > labels_;

	std::string build_plot_script_common(
			std::string quantityLabel,
			std::pair<double,double> minMax) const;

	std::string build_xtics() const;

	void dump_file(
			std::string const & filename,
			std::string const & content) const;

};

} /* namespace IOMethods */
} /* namespace elephon */

#endif /* ELEPHON_IOMETHODS_KPATH_H_ */
