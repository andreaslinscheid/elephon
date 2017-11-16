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
#include <memory>

namespace elephon
{
namespace PhononStructure
{

class Phonon
{
public:

	void initialize(std::shared_ptr<const ForceConstantMatrix> fc,
			std::vector<double> masses);

	/**
	 * Diagonalize the Fourier transformed matrix of force constants at a list of q vectors.
	 *
	 * @param q				List of q vectors in the layout q1x, q1y, q1z, q2x ...
	 * @param w2			A vector that will be resized to fit the phonon frequencies in units of THz.
	 * 						Layout is q1 mode 0, mod 1 ... mode nM, q2 ...
	 * @param eigenModes	A vector that will be resized to fit the dynamical matrices.
	 * 						Layout is (q, mode mu, mode nu) such that nu is fastest and q is slowest running.
	 */
	void compute_at_q(std::vector<double> const & q,
			std::vector<double> & w2,
			std::vector< std::complex<double> > & eigenModes) const;

	void evaluate_derivative(
			std::vector<double> const & q,
			std::vector<double> & dwdq) const;

	/**
	 * Other name for compute_at_q to conform with templated band structure calculations.
	 */
	void evaluate(std::vector<double> const & q,
			std::vector<double> & w2,
			std::vector< std::complex<double> > & eigenModes) const;

	int get_num_modes() const;

	std::vector<double> const & get_masses() const;

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
