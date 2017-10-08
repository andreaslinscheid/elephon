/*	This file LocalDerivatives.hpp is part of elephon.
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
 *  Created on: Sep 28, 2017
 *      Author: A. Linscheid
 */

#include "Algorithms/LocalDerivatives.h"
#include "Algorithms/LinearAlgebraInterface.h"

namespace elephon
{
namespace Algorithms
{
namespace localDerivatives
{

template<typename T, class DataLoader>
void
compute_derivatives_sqr_polynom(
		int nBnd,
		std::vector<int> const & reducibleKPTIndices,
		std::vector<T> * gradientFieldPtr,
		std::vector<T> * hessianFieldPtr,
		LatticeStructure::RegularBareGrid const & grid,
		DataLoader const & reducibleData)
{
	if ( (nBnd == 0)
			or reducibleKPTIndices.empty()
			or ((gradientFieldPtr == nullptr) and (hessianFieldPtr == nullptr)) )
		return;

	auto const & d = grid.get_grid_dim();
	elephon::Algorithms::LinearAlgebraInterface linalg;

	// construct the cubic model pseudo inverse matrix for least square fit
	int nKn = 19;
	std::vector<double> k(3);
	int nCoeff = 19;
	std::vector<double> A(nCoeff*nKn, 1);
	int n = 0;
	for ( int ipz = -1; ipz <= 1; ++ipz)
		for ( int ipy = -1; ipy <= 1; ++ipy)
			for ( int ipx = -1; ipx <= 1; ++ipx)
			{
				if ( abs(ipx)+abs(ipy)+abs(ipz) == 3 )
					continue;
				k[0] = ipx/double(d[0]);
				k[1] = ipy/double(d[1]);
				k[2] = ipz/double(d[2]);
				// go to cartesian coords in units of inverse angstroems
				grid.get_lattice().reci_direct_to_cartesian_2pibya(k);

				double x = k[0];
				double y = k[1];
				double z = k[2];
				int c = 0;
				// constant excluded
				A[nCoeff*n + c++] = x;
				A[nCoeff*n + c++] = y;
				A[nCoeff*n + c++] = z;
				A[nCoeff*n + c++] = x*x;
				A[nCoeff*n + c++] = 2*x*y;
				A[nCoeff*n + c++] = 2*x*z;
				A[nCoeff*n + c++] = y*y;
				A[nCoeff*n + c++] = 2*y*z;
				A[nCoeff*n + c++] = z*z;

				A[nCoeff*n + c++] = x*x*x;
				A[nCoeff*n + c++] = 3*x*x*y;
				A[nCoeff*n + c++] = 3*x*x*z;
				A[nCoeff*n + c++] = 3*x*y*y;
				A[nCoeff*n + c++] = 6*x*y*z;
				A[nCoeff*n + c++] = 3*x*z*z;

				A[nCoeff*n + c++] = y*y*y;
				A[nCoeff*n + c++] = 3*y*y*z;
				A[nCoeff*n + c++] = 3*y*z*z;

				A[nCoeff*n + c++] = z*z*z;
				n++;
			}
	std::vector<double> AInv;
	linalg.pseudo_inverse(std::move(A), nKn, nCoeff, AInv );

	std::vector<double> bval(nKn);
	std::vector<double> fit(nCoeff);

	if ( gradientFieldPtr  != nullptr)
		gradientFieldPtr->resize(nBnd*reducibleKPTIndices.size()*3);
	if ( hessianFieldPtr  != nullptr)
		hessianFieldPtr->resize(nBnd*reducibleKPTIndices.size()*6);

	for ( int ib = 0 ; ib < nBnd; ++ib)
		for ( int ikrm = 0 ; ikrm < reducibleKPTIndices.size(); ++ikrm)
		{
			int ikr = reducibleKPTIndices[ikrm];
			auto xyz = grid.get_reducible_to_xyz(ikr);
			auto xyzMod = xyz;
			int n = 0;
			for ( int ipz = -1; ipz <= 1; ++ipz)
				for ( int ipy = -1; ipy <= 1; ++ipy)
					for ( int ipx = -1; ipx <= 1; ++ipx)
					{
						if ( abs(ipx)+abs(ipy)+abs(ipz) == 3 )
							continue;
						xyzMod[0] = xyz[0] + ipx;
						xyzMod[1] = xyz[1] + ipy;
						xyzMod[2] = xyz[2] + ipz;

						for ( int i = 0 ; i < 3 ; ++i)
						{
							xyzMod[i] = xyzMod[i] < 0 ? xyzMod[i] + d[i] : xyzMod[i];
							xyzMod[i] = xyzMod[i] >= d[i] ? xyzMod[i] - d[i] : xyzMod[i];
						}
						int ikrn = grid.get_xyz_to_reducible(xyzMod);
						bval[n] = reducibleData(ikrn, ib)-reducibleData(ikr, ib);
						n++;
					}

			std::fill(fit.begin(), fit.end(), 0.0);
			for ( int ip = 0 ; ip < nCoeff; ++ip)
				for ( int ikn = 0 ; ikn < nKn; ++ikn)
					fit[ip] += AInv[ip*nKn+ikn]*bval[ikn];

			int offset = ikrm*nBnd+ib;
			if ( gradientFieldPtr  != nullptr)
			{
				auto it = gradientFieldPtr->begin() + offset*3;
				std::copy(&fit[0], &fit[0]+3, it);
			}

			// note: the factor 2 is from d^2 (x^2)/dx^2 = 2
			if ( hessianFieldPtr  != nullptr)
			{
				auto it = hessianFieldPtr->begin() + offset*6;
				*(it+0) = 2*fit[3+0]; // x x
				*(it+1) =   fit[3+1]; // x y
				*(it+2) =   fit[3+2]; // x z
				*(it+3) = 2*fit[3+3]; // y y
				*(it+4) =   fit[3+4]; // y z
				*(it+5) = 2*fit[3+5]; // z z
			}
		}
}

} /* namespace localDerivatives */
} /* namespace Algorithms */
} /* namespace elephon */
