/*	This file Phonon.cpp is part of elephon.
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

#include "PhononStructure/Phonon.h"
#include "Algorithms/LinearAlgebraInterface.h"

namespace elephon
{
namespace PhononStructure
{

void
Phonon::initialize( ForceConstantMatrix fc,
		std::vector<double> masses)
{
	fc_ = std::move(fc);
	masses_ = std::move(masses);
	assert( fc_.get_num_modes() == 3*masses_.size() );
}

void
Phonon::compute_at_q(std::vector<double> const & q,
		std::vector<double> & w2,
		std::vector< std::complex<double> > & eigenModes) const
{
	assert( q.size()%3 == 0 );
	int nq = q.size()/3;
	int nM = fc_.get_num_modes();

	Algorithms::LinearAlgebraInterface linalg;
	std::vector< std::complex<double> > qlocalModes( nM*nM ),
			 qlocalBuf; //used for swapping
	std::vector<double> qlocalFreq;

	//construct a mass array - we use that the mu layout is (atom,[x,y,z])
	assert(masses_.size() == nM/3);
	std::vector<double> sqrtMasses(nM);
	for ( int mu1 = 0 ; mu1 < nM ; ++mu1)
		sqrtMasses[mu1] = std::sqrt(masses_[mu1/3]);

	fc_.fourier_transform_q(q,eigenModes);
	assert( eigenModes.size() == nM*nM*nq );
	w2.resize( nq*nM );
	for ( int iq = 0 ; iq < nq ; ++iq)
	{
		for ( int mu1 = 0 ; mu1 < nM ; ++mu1)
			for ( int mu2 = 0 ; mu2 < nM ; ++mu2)
				qlocalModes[mu1*nM+mu2] = eigenModes[ (iq*nM+mu1)*nM+mu2 ]
						 / sqrtMasses[mu1] / sqrtMasses[mu2] ;

		linalg.diagonalize_hermitian( true, true, std::move(qlocalModes), qlocalBuf, qlocalFreq);
		std::swap(qlocalBuf,qlocalModes);

		std::copy( qlocalModes.begin(), qlocalModes.end(), eigenModes.begin() + iq*nM*nM );
		for ( int mu = 0 ; mu < nM ; ++mu)
			w2[iq*nM+mu] = std::sqrt(eVToTHzConversionFactor_)*(
					qlocalFreq[mu] >= 0 ? std::sqrt(qlocalFreq[mu]) : -std::sqrt(-qlocalFreq[mu]));
	}
}

int
Phonon::get_num_modes() const
{
	return fc_.get_num_modes();
}

std::vector<double> const &
Phonon::get_masses() const
{
	return masses_;
}

} /* namespace PhononStructure */
} /* namespace elephon */
