/*	This file Phonon.h is part of elephon.
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
 *  Created on: Jun 22, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_PHONONSTRUCTURE_PHONON_H_
#define ELEPHON_PHONONSTRUCTURE_PHONON_H_

#include "PhononStructure/ForceConstantMatrix.h"
#include "IOMethods/KPath.h"
#include "Auxillary/AlignedVector.h"
#include <memory>

namespace elephon
{
namespace PhononStructure
{

/**
 * Object to compute the Phonon modes from a matrix of force constants and atom masses.
 *
 * This object does not directly pre-compute Phonon data, but it stores information to compute
 * it (or its momentum derivative) at any given reciprocal lattice vector.
 */
class Phonon
{
public:

	/**
	 * Set up the with the data required to compute phonon dispersions.
	 *
	 * @param fc		Pointer to the force constant data.
	 * @param masses	Vector of size of the number of atoms in the primitive unit cell with their respective masses.
	 * 					Must be in the same order as the atoms occur in LatticeStructure::UnitCell
	 */
	void initialize(
			std::shared_ptr<const ForceConstantMatrix> fc,
			std::vector<double> masses);

	/**
	 * Diagonalize the Fourier transformed matrix of force constants at a list of q vectors.
	 *
	 * @param q				List of q vectors in the layout q1x, q1y, q1z, q2x ...
	 * @param w2			Resized to fit the phonon frequencies in units of THz.
	 * 						Layout is [q_index][mode_index]
	 * @param eigenModes	Resized to fit the dynamical matrices.
	 * 						Layout is [q_index][mode mu][mode nu].
	 */
	void compute_at_q(std::vector<double> const & q,
			Auxillary::Multi_array<double,2> & w2,
			Auxillary::Multi_array<std::complex<double>,3> & eigenModes) const;

	/**
	 * Similar to compute_at_q except that not the phonon mode, but its derivative.
	 *
	 * @param q			List of q vectors in the layout q1x, q1y, q1z, q2x ...
	 * @param dwdq		A vector that will be resized to fit the gradient of the phonon frequencies in units of THz.
	 */
	void evaluate_derivative(
			std::vector<double> const & q,
			Auxillary::alignedvector::DV & dwdq) const;

	/**
	 * Other name for compute_at_q to conform with templated band structure calculations.
	 */
	void evaluate(std::vector<double> const & q,
			Auxillary::Multi_array<double,2> & w2,
			Auxillary::Multi_array<std::complex<double>,3> & eigenModes) const;

	/**
	 * Get number of modes.
	 * @return	The number of modes.
	 */
	int get_num_modes() const;

	/**
	 * Get a list of the masses of atoms.
	 *
	 * @return	A vector with the masses of the atoms
	 */
	std::vector<double> const & get_masses() const;

	/**
	 * Produce a plot of the phonons along a path in reciprocal space.
	 *
	 * @param filename		Base name of the file. Will also produce a file <filename>.gp with a gnuplot script to plot it.
	 * @param kpath			The reciprocal path object.
	 */
	void write_bands_path(
			std::string const & filename,
			std::shared_ptr<const IOMethods::KPath> kpath ) const;
private:

	std::shared_ptr<const ForceConstantMatrix> fc_;

	std::vector<double> masses_;
};

} /* namespace PhononStructure */
} /* namespace elephon */

#endif /* ELEPHON_PHONONSTRUCTURE_PHONON_H_ */
