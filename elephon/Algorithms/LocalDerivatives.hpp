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
#include <set>

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
	compute_derivatives_sqr_polynom_symmetric(
			nBnd,
			reducibleKPTIndices,
			gradientFieldPtr,
			hessianFieldPtr,
			grid,
			LatticeStructure::Symmetry(),
			reducibleData);
}

template<typename T, class DataLoader>
void
compute_derivatives_sqr_polynom_symmetric(
		int nBnd,
		std::vector<int> const & reducibleKPTIndices,
		std::vector<T> * gradientFieldPtr,
		std::vector<T> * hessianFieldPtr,
		LatticeStructure::RegularBareGrid const & grid,
		LatticeStructure::Symmetry const & sym,
		DataLoader const & reducibleData)
{
	if ( (nBnd == 0)
			or reducibleKPTIndices.empty()
			or ((gradientFieldPtr == nullptr) and (hessianFieldPtr == nullptr)) )
		return;

	auto const & d = grid.get_grid_dim();

	if ( gradientFieldPtr  != nullptr)
		gradientFieldPtr->resize(nBnd*reducibleKPTIndices.size()*3);
	if ( hessianFieldPtr  != nullptr)
		hessianFieldPtr->resize(nBnd*reducibleKPTIndices.size()*6);

	for (int ikE = 0 ; ikE < reducibleKPTIndices.size(); ++ikE)
	{
		elephon::Algorithms::LinearAlgebraInterface linalg;

		std::vector<double> kE = grid.get_vector_direct(reducibleKPTIndices[ikE]);
		auto smallGroupK = sym;
		smallGroupK.set_reciprocal_space_sym(true);
		auto dropSym = smallGroupK.small_group(kE);

		// generate the set of all kpoints +-1 distant from this in any x y or z
		// plus their symmetry equivalents.
		std::set<int> reducibleKPts;
		auto xyz = grid.get_reducible_to_xyz(reducibleKPTIndices[ikE]);
		for ( int ipz = -1; ipz <= 1; ++ipz)
			for ( int ipy = -1; ipy <= 1; ++ipy)
				for ( int ipx = -1; ipx <= 1; ++ipx)
				{
					if ( (ipx == 0) && (ipy == 0) && (ipz == 0) )
						continue;
					auto xyzNeighbor = xyz;
					xyzNeighbor[0] += ipx;
					xyzNeighbor[1] += ipy;
					xyzNeighbor[2] += ipz;

					reducibleKPts.insert(grid.get_xyz_to_reducible_periodic(xyzNeighbor));
				}

		// expand the set with symmetry operations of the small group of kE
		std::set<int> symmExpanded;
		std::vector<int> rotatedIndices;
		for (int isymK = 0 ; isymK < smallGroupK.get_num_symmetries_no_T_rev(); ++isymK)
		{
			auto kVectors = grid.get_vectors_direct(std::vector<int>(reducibleKPts.begin(), reducibleKPts.end()));
			smallGroupK.rotate<double>(isymK, kVectors.begin(), kVectors.end(), true);
			grid.get_list_reducible_lattice_point_indices(kVectors, rotatedIndices);
			symmExpanded.insert(rotatedIndices.begin(), rotatedIndices.end());
		}

		// construct the cubic model pseudo inverse matrix for least square fit
		int nKn = symmExpanded.size();
		std::vector<double> k(3);
		int nCoeff = 19;
		std::vector<double> A(nCoeff*nKn, 1);
		int n = 0;
		std::vector<int> xyzN(3);

		std::vector<double> bval(nKn*nBnd);

		for (int ik : symmExpanded)
		{
			// take the periodic copy which is closes to the center point
			grid.get_reducible_to_xyz(ik,xyzN);
			for (int i = 0 ; i < 3; ++i)
			{
				int relDistance = xyzN[i] - xyz[i];
				if (std::abs(xyzN[i]+d[i]-xyz[i])<std::abs(relDistance))
					relDistance = xyzN[i] - xyz[i] + d[i];

				if (std::abs(xyzN[i]-d[i]-xyz[i])<std::abs(relDistance))
					relDistance = xyzN[i] - xyz[i] - d[i];

				xyzN[i] = relDistance;
				k[i] = double(xyzN[i])/double(d[0]);
			}

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

			for ( int ib = 0 ; ib < nBnd; ++ib)
				bval[ib*nKn+n] = reducibleData(ik, ib)-reducibleData(reducibleKPTIndices[ikE], ib);
			n++;
		}
		std::vector<double> AInv;
		linalg.pseudo_inverse(std::move(A), nKn, nCoeff, AInv );

		std::vector<double> fit(nCoeff);
		for ( int ib = 0 ; ib < nBnd; ++ib)
		{
			std::fill(fit.begin(), fit.end(), 0.0);
			for ( int ip = 0 ; ip < nCoeff; ++ip)
				for ( int ikn = 0 ; ikn < nKn; ++ikn)
					fit[ip] += AInv[ip*nKn+ikn]*bval[ib*nKn+ikn];

			int offset = ikE*nBnd+ib;
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
}

} /* namespace localDerivatives */
} /* namespace Algorithms */
} /* namespace elephon */
