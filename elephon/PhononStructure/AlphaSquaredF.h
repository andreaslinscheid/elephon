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
	 * This is the main driver for the a2F on the dense wavefunction grid calculation.
	 *
	 * Computes and sets the a2F coupling function internally.
	 *
	 * @param resourceHandler	The interface on how to obtain all the data need for this task.
	 */
	void compute_a2F_grid( std::shared_ptr<IOMethods::ResourceHandler> resourceHandler );

	void write_a2F_file(std::string const & filename) const;
private:

	std::vector<double> a2F_;

	double freqMin_ = 0;

	double freqMax_ = 0;

	int freqNPts_ = 0;

	void setup_internal_freq_grid(
			std::shared_ptr<const PhononStructure::PhononGrid> phgrid,
			std::vector<double> const & phrange,
			int npts);

	void map_freq_grid_slot(
			std::vector<float>::const_iterator begData,
			std::vector<float>::const_iterator begFreq,
			std::vector<float>::const_iterator endFreq);

	/**
	 * Build a list of q points on the grid where each elements lists the k and k' grid points associated with this q.
	 *
	 *  The k and k' grid points are obtained by first computing the closest grid point of the kList and kpList. Then,
	 *  these grid points are used to compute q=k'-k modulo a reciprocal lattice vector. The method computes for each q that
	 *  occurs a list of k and k' grid indices associated with this vector.
	 *
	 * @param kList				a vector with 3*nk elements 1x,1y,1z,2x ... of arbitrary (not necessarily on the grid) k points
	 * @param kpList			a vector with 3*nk elements 1x,1y,1z,2x ... of arbitrary (not necessarily on the grid) k' points
	 * @param kPointToGrid		(output) A list that will be resized to hold for each grid k grid point that occurs a pair with
	 * 									its reducible index and a list of indices in the \p kList that map to this reducible index.
	 * @param kpPointToGrid		(output) A list that will be resized to hold for each grid k' grid point that occurs a pair with
	 * 									its reducible index and a list of indices in the \p kpList that map to this reducible index
	 * @param grid				the bare regular grid where the k, k' and q points will be discretized on.
	 * @param q_to_k_and_kp_map	(output) A list that will be resized to hold for each q=k'-k point in the grid that occurs
	 * 									a pair of this q-grid point and a list of (k,k') index pairs where each index refers to an
	 * 									index in the kPointToGrid (first) and kpPointToGrid (second) list, respectively.
	 * 									The list consists of all k'-k grid points to result
	 * 									in this q vector module a reciprocal lattice vector.
	 */
	void query_q(
			std::vector<double> const & kList,
			std::vector<double> const & kpList,
			std::vector<std::pair<int,std::vector<int>>> & kPointToGrid,
			std::vector<std::pair<int,std::vector<int>>> & kpPointToGrid,
			LatticeStructure::RegularBareGrid const & grid,
			std::vector<std::pair<int, std::vector<std::pair<int,int>>>> & q_to_k_and_kp_map);
};

} /* namespace PhononStructure */
} /* namespace elephon */

#endif /* ELEPHON_PHONONSTRUCTURE_ALPHASQUAREDF_H_ */
