/*	This file ElectronPhononCouling.h is part of elephon.
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
 *  Created on: Jul 4, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_PHONONSTRUCTURE_ELECTRONPHONONCOUPLING_H_
#define ELEPHON_PHONONSTRUCTURE_ELECTRONPHONONCOUPLING_H_

#include "PhononStructure/Phonon.h"
#include "PhononStructure/DisplacementPotential.h"
#include "ElectronicStructure/Wavefunctions.h"
#include <vector>
#include <complex>

namespace elephon
{
namespace PhononStructure
{

class ElectronPhononCoupling
{
public:

	void generate_gkkp(
			std::vector<double> kList,
			std::vector<double> kpList,
			std::vector<int> bandList,
			std::vector<int> bandpList,
			Phonon const & ph,
			DisplacementPotential const & dvscf,
			ElectronicStructure::Wavefunctions const & wfcts );

	void write_gkkp_file(
			std::string const & filename,
			std::vector<double> const & kList,
			std::vector<double> const & kpList,
			std::vector<int> const & bandList,
			std::vector<int> const & bandpList) const;

	std::complex<float> operator() (int ik, int ikp, int ib, int ibp, int imu) const;

	void get_local_matrix_range(int ik, int ikp,
			std::vector< std::complex<float> >::iterator & rangeBegin,
			std::vector< std::complex<float> >::iterator & rangeEnd );

private:

	int nK_ = 0;

	int nB_ = 0;

	int nKp_ = 0;

	int nBp_ = 0;

	int nM_ = 0;

	std::vector<std::complex<float>> data_;

	int tensor_layout(int ik, int ikp, int ib, int ibp, int imu) const;
};

} /* namespace PhononStructure */
} /* namespace elephon */

#endif /* ELEPHON_PHONONSTRUCTURE_ELECTRONPHONONCOUPLING_H_ */
