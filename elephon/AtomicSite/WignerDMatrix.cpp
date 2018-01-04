/*	This file WignerDMatrix.cpp is part of elephon.
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
 *  Created on: Jan 2, 2018
 *      Author: A. Linscheid
 */

#include "AtomicSite/WignerDMatrix.h"
#include "Algorithms/helperfunctions.hpp"

namespace elephon
{
namespace AtomicSite
{

void
WignerDMatrix::initialize(
		int angularQuantumNumber,
		double alpha,
		double beta,
		double gamma)
{
	l_ = angularQuantumNumber;
	dim_ = 2*l_+1;

	// map beta to range [0,2pi[
	beta -= 2.0*M_PI*std::floor(beta/(2.0*M_PI));

	// initialize with the real small d matrix
	std::vector<double> smallDMatrix;
	this->wigner_small_d(beta, smallDMatrix);
	assert(smallDMatrix.size() == dim_*dim_);
	matrix_.assign(smallDMatrix.begin(), smallDMatrix.end());

	for ( int m = -l_; m <= l_; ++m )
		for ( int mp = -l_; mp <= l_; ++mp )
		{
			matrix_[this->angular_layout(m,mp)] *= std::complex<double>(std::cos(alpha*m), std::sin(alpha*m))
							*std::complex<double>(std::cos(gamma*mp), std::sin(gamma*mp));
		}
}

void
WignerDMatrix::wigner_small_d(
		double beta,
		std::vector<double> & D) const
{
	// https://doi.org/10.1063/1.2194548 (Eq. 30)
	if ( (beta > 0) && (beta <= M_PI/2.0) )
	{
		this->wigner_small_d_first_quad(beta, D);
	}
	else if ( beta == 0 )
	{
		// delta m, k
		D.assign(dim_*dim_, 0.0);
		for (int i = 0 ; i < dim_; ++i )
			D[i*dim_+i] = 1.0;
	}
	else if ( (beta > M_PI/2.0) && (beta < M_PI))
	{
		this->wigner_small_d_first_quad(M_PI-beta, D);
		auto Dp = D;
		for (int m = -l_; m <= l_; ++m)
			for (int k = -l_; k <= l_; ++k)
				D[this->angular_layout(m,k)] = this->power_minus_1(l_+k)*Dp[this->angular_layout(-m,k)];
	}
	else if ( beta == M_PI)
	{
		D.assign(dim_*dim_, 0.0);
		for (int m = -l_; m <= l_; ++m)
			for (int k = -l_; k <= l_; ++k)
				D[this->angular_layout(m,k)] = this->power_minus_1(l_+k)*(-m == k ? 1.0 : 0.0);
	}
	else if ( (beta > M_PI) && (beta < 3.0*M_PI/2.0))
	{
		this->wigner_small_d_first_quad(beta-M_PI, D);
		auto Dp = D;
		for (int m = -l_; m <= l_; ++m)
			for (int k = -l_; k <= l_; ++k)
				D[this->angular_layout(m,k)] = this->power_minus_1(l_+k)*Dp[this->angular_layout(m,-k)];
	}
	else if ( (beta >= 3.0*M_PI/2.0) && (beta < 2.0*M_PI))
	{
		this->wigner_small_d_first_quad(2.0*M_PI-beta, D);
		for (int m = -l_; m <= l_; ++m)
			for (int k = -l_; k <= l_; ++k)
				D[this->angular_layout(m,k)] *= this->power_minus_1(m+k);
	}
	else
	{
		throw std::logic_error("Called wigner_small_d with incorrect angle");
	}
}

void
WignerDMatrix::wigner_small_d_first_quad(
		double beta,
		std::vector<double> & D) const
{
	assert( (beta > 0) and (beta <= M_PI/2.0) );
	D.resize(dim_*dim_);

	double sin_b = std::sin(beta);
	double cos_b = std::cos(beta);
	auto glm_mat_order = [&] (int m, int mp) { return m*(l_+1)+mp;};

	// https://doi.org/10.1063/1.2194548 (Eq. 29)
	std::vector<double> glm(std::pow(l_+1, 2), 1.0);
	for(int m = 1; m <= l_; ++m)
	{
		glm[glm_mat_order(m,0)] = std::sqrt(static_cast<double>(2*m-1)/(2*m)) * glm[glm_mat_order(m-1,0)];
		for(int mp = 1; mp <= l_; ++mp)
			glm[glm_mat_order(m,mp)] = std::sqrt(static_cast<double>(m-mp+1)/(m+mp))*glm[glm_mat_order(m,mp-1)];
	}

	// https://doi.org/10.1063/1.2194548 (Eq. 28)
	for(int m=0; m <= l_; ++m)
		D[this->angular_layout(m,l_)] = this->power_minus_1(l_+m)*glm[glm_mat_order(l_,m)]
									*std::pow(1.0+cos_b,m) * std::pow(sin_b,l_-m);

	// https://doi.org/10.1063/1.2194548 (Eq. 26)
	for(int k = l_; k > -l_; --k)
		D[this->angular_layout(l_,k-1)] = (l_+k) / std::sqrt(static_cast<double>(l_*(l_+1)-k*(k-1)))
									*(sin_b/(1.0+cos_b))* D[this->angular_layout(l_,k)];

	// https://doi.org/10.1063/1.2194548 (Eq. 25)
	for(int m = l_-1; m >= 0; --m)
	  for(int k = l_; k > -l_; --k)
		  D[this->angular_layout(m,k-1)] = std::sqrt(static_cast<double>(l_*(l_+1)-m*(m+1)))
									  / std::sqrt(static_cast<double>(l_*(l_+1)-k*(k-1))) * D[this->angular_layout(m+1,k)]
									+ (m+k) / std::sqrt(static_cast<double>(l_*(l_+1)-k*(k-1)))
										*sin_b/(1.0+cos_b) * D[this->angular_layout(m,k)];

	// Symmetry (Eq. 27) to obtain negative m
	for (int m = -l_; m < 0; ++m)
	  for (int k = -l_; k <= l_; ++k)
		  D[this->angular_layout(m,k)] = this->power_minus_1(m+k) * D[this->angular_layout(-m,-k)];
}


int
WignerDMatrix::angular_layout(int m, int mp) const
{
	assert((m>=-l_)&&(m<=l_));
	assert((mp>=-l_)&&(mp<=l_));
	int index = (l_+m)*dim_+(l_+mp);
	assert(index >= 0);
	assert(dim_*dim_>index);
	return index;
}

std::complex<double>
WignerDMatrix::operator() (int m, int mp) const
{
	return matrix_[this->angular_layout(m,mp)];
}

double
WignerDMatrix::power_minus_1 (int m) const
{
	return m%2 == 0 ? 1.0 : -1.0;
};

} /* namespace AtomicSite */
} /* namespace elephon */
