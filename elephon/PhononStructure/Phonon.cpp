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
#include "Auxillary/UnitConversion.h"
#include <stdexcept>
#include <utility>

namespace elephon
{
namespace PhononStructure
{

void
Phonon::initialize( std::shared_ptr<const ForceConstantMatrix> fc,
		std::vector<double> masses)
{
	fc_ = fc;
	masses_ = std::move(masses);
	assert( fc_->get_num_modes() == 3*masses_.size() );
}

void
Phonon::compute_at_q(std::vector<double> const & q,
		Auxillary::alignedvector::DV & w2,
		Auxillary::alignedvector::ZV & eigenModes) const
{
	assert( fc_ );
	assert( q.size()%3 == 0 );
	int nq = q.size()/3;
	int nM = fc_->get_num_modes();

	Algorithms::LinearAlgebraInterface linalg;
	std::vector< std::complex<double> > qlocalModes( nM*nM ),
			 qlocalBuf; //used for swapping
	std::vector<double> qlocalFreq;

	//construct a mass array - we use that the mu layout is (atom,[x,y,z])
	assert(masses_.size() == nM/3);
	std::vector<double> sqrtMasses(nM);
	for ( int mu1 = 0 ; mu1 < nM ; ++mu1)
		sqrtMasses[mu1] = std::sqrt(masses_[mu1/3]);

	fc_->fourier_transform_q(q, eigenModes, true);
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
			w2[iq*nM+mu] = Auxillary::units::SQRT_EV_BY_A2_U_TO_THZ*
					( qlocalFreq[mu] >= 0 ? std::sqrt(qlocalFreq[mu]) : -std::sqrt(-qlocalFreq[mu]));

		// phase convention: choose the first eigenmode real
		for (int imode = 0 ; imode < nM; ++imode)
		{
			std::complex<double> phase(1.0);
			for (int icomponent = 0 ; icomponent < nM; ++icomponent)
			{
				if ( std::abs(eigenModes[(iq*nM+icomponent)*nM+imode]) > 1e-8 )
				{
					phase = eigenModes[(iq*nM+icomponent)*nM+imode]
									   /std::abs(eigenModes[(iq*nM+icomponent)*nM+imode]);
					break;
				}
			}
			for (int icomponent = 0 ; icomponent < nM; ++icomponent)
				eigenModes[(iq*nM+icomponent)*nM+imode] *= std::conj(phase);
		}
	}
}

void
Phonon::evaluate_derivative(
		std::vector<double> const & q,
		Auxillary::alignedvector::DV & dwdq) const
{
	assert( fc_ );
	Auxillary::alignedvector::DV w;
	Auxillary::alignedvector::ZV unitaryTrafo;
	this->compute_at_q(q, w, unitaryTrafo);
	Auxillary::alignedvector::ZV ftderivative;
	fc_->fourier_transform_derivative(q,ftderivative);

	int nq = q.size()/3;
	int nM = fc_->get_num_modes();

	assert(masses_.size() == nM/3);
	std::vector<double> sqrtMasses(nM);
	for ( int mu1 = 0 ; mu1 < nM ; ++mu1)
		sqrtMasses[mu1] = std::sqrt(masses_[mu1/3]);

	for ( int iq = 0 ; iq < nq ; ++iq)
		for ( int mu1 = 0 ; mu1 < nM ; ++mu1)
			for ( int mu2 = 0 ; mu2 < nM ; ++mu2)
				for ( int i = 0 ; i < 3 ; ++i)
					ftderivative[((iq*3+i)*nM+mu1)*nM+mu2] *=
							std::pow(Auxillary::units::SQRT_EV_BY_A2_U_TO_THZ,2)/(sqrtMasses[mu1] * sqrtMasses[mu2]);

	Algorithms::LinearAlgebraInterface linalg;

	Auxillary::alignedvector::ZV matrixBuffer(nM*nM);
	std::vector< std::complex<double> > matrixBufferResult(nM*nM);
	dwdq.resize( 3*nq*nM );
	for ( int iq = 0 ; iq < nq ; ++iq)
	{
		for ( int i = 0 ; i < 3 ; ++i )
		{
			// M = (U*) . C
			linalg.call_gemm(
					'c','n',
					nM, nM, nM, std::complex<double>(1.0),
					&unitaryTrafo[iq*nM*nM],
					nM,
					&ftderivative[(iq*3+i)*nM*nM], nM,
					std::complex<double>(0.0),
					matrixBuffer.data(),nM);

			// R = M . U
			linalg.call_gemm(
					'n','n',
					nM, nM, nM, std::complex<double>(1.0),
					matrixBuffer.data(),
					nM,
					&unitaryTrafo[iq*nM*nM], nM,
					std::complex<double>(0.0),
					matrixBufferResult.data(),nM);

			double imagSum = 0;
			for ( int inu = 0 ; inu < nM ; ++inu )
			{
				dwdq[(iq*nM+inu)*3+i] = std::real(matrixBufferResult[inu*nM+inu]) / (2.0 * w[iq*nM+inu]);
				imagSum += std::imag(matrixBufferResult[inu*nM+inu]);
			}
			if ( imagSum > 1e-4 )
				throw std::logic_error(std::string("Imaginary part of derivative is ")
										+std::to_string(imagSum)+ ", too large");
		}
	}
}

void
Phonon::evaluate(std::vector<double> const & q,
		Auxillary::alignedvector::DV & w2,
		Auxillary::alignedvector::ZV & eigenModes) const
{
	this->compute_at_q(q, w2, eigenModes);
}

int
Phonon::get_num_modes() const
{
	return fc_->get_num_modes();
}

std::vector<double> const &
Phonon::get_masses() const
{
	return masses_;
}

void
Phonon::write_bands_path(
		std::string const & filename,
		std::shared_ptr<const IOMethods::KPath> kpath ) const
{
	auto gnuplotFile = filename+".gp";

	Auxillary::alignedvector::DV modesAlongPath;
	std::vector<double> const & qpath = kpath->get_k_points();
	Auxillary::alignedvector::ZV dynmat;
	this->compute_at_q(qpath, modesAlongPath, dynmat);

	auto mm = std::minmax_element(modesAlongPath.begin(), modesAlongPath.end());
	double range = *mm.second - *mm.first;
	auto frequencyWindow = std::make_pair(*mm.first, *mm.second+range*0.05);

	kpath->produce_gnuplot_script_stable_particle(
			gnuplotFile,
			filename,
			"${\varOmega}_{\nu}(\bf{q})$",
			modesAlongPath,
			this->get_num_modes(),
			frequencyWindow);
}

} /* namespace PhononStructure */
} /* namespace elephon */
