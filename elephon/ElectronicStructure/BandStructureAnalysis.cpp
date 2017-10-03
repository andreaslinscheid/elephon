/*	This file BandStructureAnalysis.cpp is part of elephon.
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
 *  Created on: Sep 7, 2017
 *      Author: A. Linscheid
 */

#include "ElectronicStructure/BandStructureAnalysis.h"
#include "Algorithms/LinearAlgebraInterface.h"
#include "Algorithms/FFTInterface.h"
#include <fstream>
#include <set>
#include <stdexcept>

#include "ElectronicStructure/GradientFFTReciprocalGrid.h"

namespace elephon
{
namespace ElectronicStructure
{
namespace BandStructureAnalysis
{

void find_band_extrema(
		ElectronicBands const & bands,
		std::vector<b_extrema> & extrema,
		std::vector<double> const & energyWindow)
{
	auto const & g = bands.get_grid();
	extrema.clear();

	auto get_neighbours = [] (
			LatticeStructure::RegularSymmetricGrid const & grid,
			int ik)
	{
		std::vector<int> neigbours(26);
		int ikr = grid.get_maps_irreducible_to_reducible()[ik][0];
		auto xyz = grid.get_reducible_to_xyz(ikr);
		auto xyzMod = xyz;
		int c = 0;
		for ( int pz = -1; pz <= 1; ++pz)
			for ( int py = -1; py <= 1; ++py)
				for ( int px = -1; px <= 1; ++px)
				{
					if ( (px == 0) and (py == 0) and (pz == 0) )
						continue;
					xyzMod[0] = xyz[0]+px;
					xyzMod[1] = xyz[1]+py;
					xyzMod[2] = xyz[2]+pz;
					for ( int i = 0 ; i < 3 ; ++i )
					{
						xyzMod[i] = xyzMod[i] < grid.get_grid_dim()[i] ? xyzMod[i] : xyzMod[i] - grid.get_grid_dim()[i];
						xyzMod[i] = xyzMod[i] >= 0 ? xyzMod[i] : xyzMod[i] + grid.get_grid_dim()[i];
					}
					neigbours[c++] = grid.get_maps_red_to_irreducible()[grid.get_xyz_to_reducible(xyzMod)];
				}
		return neigbours;
	};

	auto indicesInWindow = bands.get_bands_crossing_energy_window(energyWindow);
	std::vector<double> compareValues(27);
	for ( int ik = 0 ; ik < g.get_np_irred(); ++ik)
	{
		auto neighbours = get_neighbours(g, ik);
		for ( int ib : indicesInWindow)
		{
			compareValues[0] = bands(ik, ib);
			if ( (compareValues[0] < energyWindow[0]) or (compareValues[0] > energyWindow[1]) )
				continue;

			for ( int in = 0 ; in < neighbours.size(); ++in )
				compareValues[1+in] = bands(neighbours[in], ib);

			auto m = std::minmax_element(compareValues.begin(), compareValues.end());
			if ( std::distance(m.first, compareValues.begin()) == 0 )
			{
				b_extrema e;
				e.ibnd = ib;
				e.k = g.get_vector_direct(g.get_maps_irreducible_to_reducible()[ik][0]);
				e.energy = compareValues[0];
				extrema.push_back(std::move(e));
			}

			// the code for maxima is a negative band index
			if ( std::distance(m.second, compareValues.begin()) == 0 )
			{
				b_extrema e;
				e.ibnd = -ib-1;
				e.k = g.get_vector_direct(g.get_maps_irreducible_to_reducible()[ik][0]);
				e.energy = compareValues[0];
				extrema.push_back(std::move(e));
			}
		}
	}
}

void compute_mass_tensor_at_extrema_poly(
		ElectronicBands const & bands,
		std::vector<b_extrema> const & extrema,
		std::vector<std::vector<double>> & massTensor,
		std::vector<double> const & energyWindow)
{
	// get a list of k points and band indices
	std::vector<double> kvects(3*extrema.size());
	std::vector<int> bandIndicesExtrema(extrema.size());
	for ( int ie = 0 ; ie < extrema.size(); ++ie)
	{
		for ( int i = 0 ; i < 3; ++i)
			kvects[ie*3 + i] = extrema[ie].k[i];
		bandIndicesExtrema[ie] = extrema[ie].ibnd < 0 ? -1 - extrema[ie].ibnd : extrema[ie].ibnd;
	}
	std::vector<int> kIndicesExtrema;
	auto const & kgrid = bands.get_grid();
	kgrid.get_list_reducible_lattice_point_indices(kvects, kIndicesExtrema);
	auto d = kgrid.get_grid_dim();

	std::set<int> bandSet(bandIndicesExtrema.begin(), bandIndicesExtrema.end());
	std::vector<float> hessian;
	bands.compute_derivatives_sqr_polynom<float>(
			std::vector<int>(bandSet.begin(), bandSet.end()),
			kIndicesExtrema,
			nullptr,
			&hessian );
	assert( hessian.size() == 6*bandSet.size()*kIndicesExtrema.size() );

	std::set<int> kptSet(kIndicesExtrema.begin(), kIndicesExtrema.end());
	massTensor.resize( extrema.size() );
	for ( int ie = 0 ; ie < extrema.size(); ++ie)
	{
		int ikr = kIndicesExtrema[ie];
		int ib = bandIndicesExtrema[ie];
		int ibm = std::distance(bandSet.begin(), bandSet.find(ib));
		int ikm = std::distance(kptSet.begin(), kptSet.find(ikr));
		assert(ibm < bandSet.size());
		assert(ikm < kptSet.size());

		// note: the factor 2 is from d^2 (x^2)/dx^2 = 2
		massTensor[ie].resize(9);
		int offset = (ikm*bandSet.size()+ibm)*6;
		massTensor[ie][0*3+0] = hessian[offset+0];
		massTensor[ie][1*3+0] = massTensor[ie][0*3+1] = hessian[offset+1];
		massTensor[ie][2*3+0] = massTensor[ie][0*3+2] = hessian[offset+2];
		massTensor[ie][1*3+1] = hessian[offset+3];
		massTensor[ie][2*3+1] = massTensor[ie][1*3+2] = hessian[offset+4];
		massTensor[ie][2*3+2] = hessian[offset+5];
	}
}


void compute_mass_tensor_at_extrema_fft(
		ElectronicBands const & bands,
		std::vector<b_extrema> const & extrema,
		std::vector<std::vector<double>> & massTensor,
		std::vector<double> const & energyWindow)
{
	auto const & g  =  bands.get_grid().get_grid_dim();
	auto indicesInWindow = bands.get_bands_crossing_energy_window(energyWindow);
	std::vector<double> bndData, massTensFullGrid;
	bands.generate_reducible_data(indicesInWindow, bndData);

	auto A = bands.get_grid().get_lattice().get_latticeMatrix();
	for ( auto &a : A)
		a *= bands.get_grid().get_lattice().get_alat();
	Algorithms::FFTInterface fft;
	fft.fft_hessian(
			A,
			g,
			bndData,
			massTensFullGrid,
			indicesInWindow.size());

	// get a list of k points and band indices
	std::vector<double> kvects(3*extrema.size());
	std::vector<int> bandIndicesExtrema(extrema.size());
	for ( int ie = 0 ; ie < extrema.size(); ++ie)
	{
		for ( int i = 0 ; i < 3; ++i)
			kvects[ie*3 + i] = extrema[ie].k[i];
		bandIndicesExtrema[ie] = extrema[ie].ibnd < 0 ? -1 - extrema[ie].ibnd : extrema[ie].ibnd;
	}

	std::vector<int> kIndicesExtrema;
	bands.get_grid().get_list_reducible_lattice_point_indices(kvects, kIndicesExtrema);

	massTensor.resize( extrema.size() );
	for ( int ie = 0 ; ie < extrema.size(); ++ie)
	{
		int cnsq = 9*(bandIndicesExtrema[ie]+indicesInWindow.size()*kIndicesExtrema[ie]);
		assert( (cnsq >=0) && (cnsq < massTensFullGrid.size()) );
		massTensor[ie].resize(9);
		std::copy(&massTensFullGrid[cnsq], &massTensFullGrid[cnsq]+9, massTensor[ie].data());
	}
}

void write_mass_tensor_file(
		std::string const & filename,
		ElectronicBands const & bands,
		int startBand,
		std::vector<double> const & energyWindow,
		int method)
{
	std::ofstream file(filename.c_str(), std::ios::binary);
	if ( ! file.good() )
		throw std::runtime_error(std::string("Problem opening file ") + filename + " for writing.");

	std::vector<b_extrema> extrema;
	find_band_extrema(bands, extrema, energyWindow);

	std::vector<std::vector<double>> massTensor;
	if ( method == 0 )
		compute_mass_tensor_at_extrema_poly(bands, extrema, massTensor, energyWindow);
	else
		compute_mass_tensor_at_extrema_fft(bands, extrema, massTensor, energyWindow);

	elephon::Algorithms::LinearAlgebraInterface lin;

	int nElem = 9+3+3+2+1;
	std::vector<float> data(nElem*extrema.size());
	std::vector<double> ev(3);
	for ( int i = 0 ; i < extrema.size(); ++i )
	{
		// Layout: first is the flag for minima (<0) or maxima (>=0)
		data[i*nElem+0] = ( extrema[i].ibnd < 0 ? 1 : -1 ) ;
		// then comes the band
		if ( extrema[i].ibnd < 0 )
			data[i*nElem+1] = startBand - 1 - extrema[i].ibnd;
		else
			data[i*nElem+1] = startBand + extrema[i].ibnd;
		// then comes the k point
		assert(extrema[i].k.size() == 3);
		data[i*nElem+2] = extrema[i].k[0];
		data[i*nElem+3] = extrema[i].k[1];
		data[i*nElem+4] = extrema[i].k[2];

		// now the energy
		data[i*nElem+5] = extrema[i].energy;

		lin.call_syev('V', 'U', 3, massTensor[i].data(), 3, ev.data() );
		// we start with the smallest eigenvalue / vector
		data[i*nElem+6 ] = ev[0];
		data[i*nElem+7 ] = massTensor[i][0*3+0];
		data[i*nElem+8 ] = massTensor[i][1*3+0];
		data[i*nElem+9 ] = massTensor[i][2*3+0];

		// middle eigenvalue / vector
		data[i*nElem+10] = ev[1];
		data[i*nElem+11] = massTensor[i][0*3+1];
		data[i*nElem+12] = massTensor[i][1*3+1];
		data[i*nElem+13] = massTensor[i][2*3+1];

		// largest eigenvalue / vector
		data[i*nElem+14] = ev[2];
		data[i*nElem+15] = massTensor[i][0*3+2];
		data[i*nElem+16] = massTensor[i][1*3+2];
		data[i*nElem+17] = massTensor[i][2*3+2];
	}

	file.write( reinterpret_cast<char*>(data.data()), sizeof(float)*data.size());
	file.close();
}

void do_band_structure_analysis(std::shared_ptr<IOMethods::ElectronicStructureCodeInterface> loader)
{
	auto const & opt = loader->get_optns();
	// Find the locations of extrema and print the mass tensor if desired ...
	if ( not opt.get_f_mtens().empty() )
	{
		ElectronicStructure::ElectronicBands bands;
		loader->read_band_structure( opt.get_root_dir(), bands);

		std::vector<double> eneWin = loader->get_optns().get_ewinbnd();
		if ( eneWin.empty() )
		{
			auto mm = bands.get_min_max();
			eneWin = std::vector<double>{mm.first, mm.second};
		}
		auto bandIndicesInWindow = bands.get_bands_crossing_energy_window( eneWin );
		auto mm = std::minmax_element(bandIndicesInWindow.begin(), bandIndicesInWindow.end() );
		auto subsetBands = bands.fft_interpolate_part(*(mm.first), *(mm.second), opt.get_fftd(), opt.get_ffts());

		int dmeth = 0;
		if ( opt.get_dmeth().compare("fft") == 0 )
			dmeth = 1;
		else if ( opt.get_dmeth().compare("pol") != 0 )
			throw std::runtime_error("unrecognized derivative method");
		write_mass_tensor_file(opt.get_f_mtens(), subsetBands, *(mm.first), eneWin, dmeth);
	}
}

} /* namespace BandStructureAnalysis */
} /* namespace ElectronicStructure */
} /* namespace elephon */
