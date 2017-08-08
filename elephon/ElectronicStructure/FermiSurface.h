/*	This file FermiSurface.h is part of elephon.
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
 *  Created on: Apr 26, 2017
 *      Author: A. Linscheid
 */

#ifndef ELECTRONICSTRUCTURE_FERMISURFACE_H_
#define ELECTRONICSTRUCTURE_FERMISURFACE_H_

#include "LatticeStructure/LatticeModule.h"
#include <vector>
#include <cstdlib>

namespace elephon
{
namespace ElectronicStructure
{

class FermiSurface
{
public:
	FermiSurface();

	void triangulate_surface(
			std::vector<int> kgrid,
			LatticeStructure::LatticeModule const& lattice,
			int nbnd,
			std::vector<double> const& energies,
			int numTargetPoints,
			double energyVal = 0);

	int get_npts_total() const;

	void get_pt(int i, std::vector<double> & p) const;

	void get_pt_weight(int i, double & pw) const;

	/**
	 * Obtain the first index of the k point associated to this band.
	 *
	 * @param ib	The index of the band
	 * @return		The index of the x coordinate of the first k vector in the
	 * 				list of Fermi vectors associated to the band \p ib
	 */
	int get_band_offset(int ib) const;

	/**
	 * Get the 3D vectors on the isosurface.
	 *
	 * The coordinates of the vectors form the list as [x1,y1,z1,x2,y2,...]
	 * in units of the direct reciprocal lattice. Different bands are seperate parts
	 * of the list, e.g. [x1,...zN,xN+1,...zM,...] if the vectors 1-N form the surface of band
	 * 1 and the vectors N+1-M form the surface of band 2. Use \see #get_band_offset(int ib) to retrieve
	 * the index of the vector of band ib.
	 *
	 * @return const reference to the list of internally stored vectors.
	 */
	std::vector<double> const& get_Fermi_vectors() const;

	/**
	 * Get the weights for each vector on the isosurface.
	 *
	 * The weights of the vectors for the list [w1,w2...]. See \see #get_Fermi_vectors
	 * on the format of the list regarding bands. Units are to be multiplied by (2p/a)^2 to
	 * give the area in reciprocal cartesian space.
	 *
	 * @return const reference to the list of internally stored weights
	 */
	std::vector<double> const& get_Fermi_weights() const;

	std::vector<double> get_Fermi_vectors_for_band(int ib) const;

	std::vector<double> get_Fermi_weights_for_band(int ib) const;
private:

	std::vector<int> kgrid_;

	std::vector<double> kfPoints_;

	std::vector<double> kfWeights_;

	//contains for a band index the index i of the first k point belonging to this bands or -1
	//if there are none.
	std::vector<int> bandsMap_;

	void band_index_range(int ib, int &start, int &end) const;
};

} /* namespace ElectronicStructure */
} /* namespace elephon */

#endif /* ELECTRONICSTRUCTURE_FERMISURFACE_H_ */
