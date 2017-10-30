/*	This file LocalDensityOfStates.cpp is part of elephon.
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
 *  Created on: Jul 24, 2017
 *      Author: A. Linscheid
 */

#include "ElectronicStructure/LocalDensityOfStates.h"
#include "ElectronicStructure/FermiSurface.h"
#include "Algorithms/TrilinearInterpolation.h"
#include "LatticeStructure/UnitCell.h"
#include "Algorithms/FFTInterface.h"
#include "Algorithms/GridRotationMap.h"
#include "Algorithms/LocalDerivatives.h"
#include "IOMethods/WriteVASPRealSpaceData.h"

namespace elephon
{
namespace ElectronicStructure
{

void
LocalDensityOfStates::compute_ldos(
		std::vector<double> energies,
		std::shared_ptr<const Wavefunctions> wfcts,
		std::shared_ptr<const LatticeStructure::UnitCell> unitcell,
		int nkpointsPerSurface,
		std::shared_ptr<const LatticeStructure::RegularBareGrid> realSpaceRes,
		std::shared_ptr<const TetrahedraIsosurface> tetraIso,
		bool symmetrize)
{
	uc_ = unitcell;
	isoEnergies_ = std::move(energies);
	rsDims_ = realSpaceRes;
	Algorithms::TrilinearInterpolation triLin( wfcts->get_k_grid().view_bare_grid() );
	Algorithms::FFTInterface fft;

	std::vector< std::vector< std::complex<float> > > wfctsFs;
	std::vector< std::complex<float> > wfctsRealSpace;
	std::vector< std::vector< int > > fftMapsWfctsFs;
	std::vector<double> gradDataAtRequestedIndices;
	std::vector<int> reqestedIndices;
	std::vector<double> FermiVelocities;
	std::vector< std::complex<float> > wfctsOneBand;
	int nrs = rsDims_->get_num_points();
	ldos_.assign( isoEnergies_.size()*nrs , 0.0 );
	std::vector<double> interpolData;
	Algorithms::FFTInterface fftInt;

	for ( int ie = 0 ; ie < isoEnergies_.size(); ++ie )
	{
		for (int ibnd = 0 ; ibnd < tetraIso->get_nBnd(); ++ibnd )
		{
			std::vector<double> isoK, isoWeights;
			tetraIso->get_irreducible_iso_vector_integration_weights(ie, ibnd, isoK, isoWeights);
			if ( isoK.empty() )
				continue;

			wfcts->generate_wfcts_at_arbitray_kp(
					isoK,
					{ibnd},
					wfctsFs,
					fftMapsWfctsFs);

			for (int ikf = 0 ; ikf < isoK.size()/3; ++ikf)
			{
				//Compute the wavefunction part |psi(r)|^2
				fft.fft_sparse_data(
						fftMapsWfctsFs[ikf],
						wfcts->get_max_fft_dims(),
						wfctsFs[ikf],
						1,
						-1,
						wfctsRealSpace,
						rsDims_->get_grid_dim(),
						false, //we want C data order
						nkpointsPerSurface);

				double contrib = isoWeights[ikf]*uc_->get_lattice().get_volume()/std::pow(2*M_PI,3);
				for ( int ir = 0 ; ir < nrs; ++ir )
				{
					ldos_[ie*nrs+ir] += contrib*(std::pow(std::real(wfctsRealSpace[ir]), 2)
										       + std::pow(std::imag(wfctsRealSpace[ir]), 2));
					assert( ! std::isnan(ldos_[ie*nrs+ir]) );
				}
			}
		}
	}

	if ( symmetrize )
	{
		// create a rotation map of grid indices
		auto const & S = uc_->get_symmetry();

		int nsym = S.get_num_symmetries();
		std::vector<std::vector<int>> rotMap;
		Algorithms::compute_grid_rotation_map({0.0, 0.0, 0.0}, *rsDims_, S, rotMap);
		assert( rotMap.size() == nsym );

		for ( int ie = 0 ; ie < isoEnergies_.size(); ++ie )
		{
			std::vector<double> ldsym(nrs, 0.0);
			for ( int isym = 0 ; isym < nsym; ++isym)
			{
				assert(rotMap[isym].size() == nrs);
				for ( int ir = 0 ; ir < nrs; ++ir)
					ldsym[rotMap[isym][ir]] += ldos_[ie*nrs+ir]/double(nsym);
			}
			std::copy(ldsym.begin(), ldsym.end(), &ldos_[ie*nrs]);
		}
	}
}

void
LocalDensityOfStates::compute_ldos(
		std::vector<double> const & energies,
		std::shared_ptr<IOMethods::ResourceHandler> loader )
{
	auto uc = loader->get_primitive_unitcell_obj();
	auto interpolGrid = loader->get_interpol_reci_mesh_obj();
	auto wfcts = loader->get_wfct_obj();
	auto rsgrid = loader->get_real_space_grid_unitcell_obj();
	auto tetraIso = loader->get_tetrahedra_isosurface();

	this->compute_ldos(
			energies,
			wfcts,
			uc,
			loader->get_optns().get_numFS(),
			rsgrid,
			tetraIso,
			loader->get_optns().get_symOut());
}

void
LocalDensityOfStates::write_file( std::string const & filename, bool binary ) const
{
	IOMethods::WriteVASPRealSpaceData output;
	auto filename_mod = filename;

	int nrs = rsDims_->get_num_points();
	for ( int ie = 0 ; ie < isoEnergies_.size(); ++ie )
	{
		if ( isoEnergies_.size() > 1 )
			filename_mod = "e"+std::to_string(ie)+"_"+filename;

		output.write_file(
				filename_mod,
				std::string("LDOS at E=")+std::to_string(isoEnergies_[ie]),
				rsDims_->get_grid_dim(),
				uc_,
				std::vector<double>(&ldos_[ie*nrs],&ldos_[ie*nrs]+nrs),
				false,
				true);
	}
}

} /* namespace ElectronicStructure */
} /* namespace elephon */
