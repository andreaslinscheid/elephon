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
#include "Auxillary/UnitConversion.h"
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

std::vector<double>
read_mass_tens_file(boost::filesystem::path massTens)
{
	std::ifstream file(massTens.c_str(), std::ios::binary);
	if ( ! file.good() )
		throw std::runtime_error("cannot open mass tensor file");

	file.seekg(0, std::ios::beg);
	int size = file.tellg();
	file.seekg(0, std::ios::end);
	size = int(file.tellg()) - size;
	std::vector<char> buf(size);
	file.seekg(0, std::ios::beg);
	file.read(&buf[0], size);

	int numElem = size/sizeof(float);

	auto ptr = reinterpret_cast<float*>(buf.data());
	return std::vector<double>(ptr, ptr+numElem);
}

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
			if ( (*m.first < energyWindow[0]) or (*m.second > energyWindow[1]) )
				continue;

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
	std::vector<double> hessian;
	bands.compute_derivatives_sqr_polynom<double>(
			std::vector<int>(bandSet.begin(), bandSet.end()),
			kIndicesExtrema,
			nullptr,
			&hessian );
	assert( hessian.size() == 6*bandSet.size()*kIndicesExtrema.size() );

	Algorithms::LinearAlgebraInterface linalg;

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

		massTensor[ie].resize(9);
		int offset = (ikm*bandSet.size()+ibm)*6;
		assert(hessian.size() > (offset+5));
		std::vector<double> h(9);
		h[0*3+0] = hessian[offset+0];
		h[1*3+0] = h[0*3+1] = hessian[offset+1];
		h[2*3+0] = h[0*3+2] = hessian[offset+2];
		h[1*3+1] = hessian[offset+3];
		h[2*3+1] = h[1*3+2] = hessian[offset+4];
		h[2*3+2] = hessian[offset+5];

		linalg.inverse(std::move(h),massTensor[ie]);

		for ( auto & mij : massTensor[ie] )
			mij *= Auxillary::units::INVERSE_EV_TIMES_A2_TO_ME;
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

void
output_valenceBandMaxima_conductionBandMinima(boost::filesystem::path massTens)
{
	auto massTensorData = ElectronicStructure::BandStructureAnalysis::read_mass_tens_file(massTens);
	int numExtrema = massTensorData.size() / 18 ;
	std::vector< std::vector<double> > kPoints(numExtrema);
	std::vector< std::vector<double> > eigenvalues(numExtrema);

	// use a compare function that treats close doubles as equal
	auto cmp = [] (double a, double b) {
		if ( std::abs(a-b)>1e-3 )
			return a < b;
		return false;
		};
	std::multimap<double,int,decltype(cmp)> sortMinima(cmp);
	std::multimap<double,int,decltype(cmp)> sortMax(cmp);
	for (int ie = 0 ; ie < numExtrema ; ++ie )
	{
		kPoints[ie] = std::vector<double>(&massTensorData[18*ie+2], &massTensorData[18*ie+2]+3);
		eigenvalues[ie] = std::vector<double>{ massTensorData[18*ie+ 6],
											   massTensorData[18*ie+10],
											   massTensorData[18*ie+14] };
		double energy = massTensorData[18*ie+5];
		if ( massTensorData[18*ie+0] >= 0 )
		{
			sortMax.insert(std::make_pair(-energy, ie));
		}
		else
		{
			sortMinima.insert(std::make_pair(energy, ie));
		}
	}

	// find the valence band maxima and the conduction band minima
	std::vector<int> cbmin;
	double last = std::numeric_limits<double>::min();
	for ( auto min : sortMinima)
	{
		if ( min.first < 0 )
			continue;

		// min.first == last.first according to cmp
		if ( (!(cmp(last, min.first)) and (!cmp(min.first ,last))) or cbmin.empty() )
		{
			cbmin.push_back(min.second);
		}
		else
		{
			break;
		}
		last = min.first;
	}
	std::vector<int> vbmax;
	for ( auto max : sortMax)
	{
		if ( max.first < 0 ) // conduction bands; note negative energy convention
			continue;

		double energy = -max.first; // revert negative energy convention.

		// min.first == last.first according to cmp
		if ( (!(cmp(last, energy)) and (!cmp(energy, last))) or vbmax.empty() )
		{
			vbmax.push_back(max.second);
		}
		else
		{
			break;
		}
		last = energy;
	}

	for ( auto i : vbmax )
	{
		std::cout << "\nValence band max at k="<<
				kPoints[i][0] << ","<< kPoints[i][1] << ","<< kPoints[i][2] << " energy " << massTensorData[18*i+5]<<"\n";
		std::cout << "Mass tensor eigenvalues: "<<
				eigenvalues[i][0] << ","<< eigenvalues[i][1] << ","<< eigenvalues[i][2] <<"\n";
	}

	for ( auto i : cbmin )
	{
		std::cout << "\nConduction band min at k="<<
				kPoints[i][0] << ","<< kPoints[i][1] << ","<< kPoints[i][2] << " energy " << massTensorData[18*i+5]<<"\n";
		std::cout << "Mass tensor eigenvalues: "<<
				eigenvalues[i][0] << ","<< eigenvalues[i][1] << ","<< eigenvalues[i][2] <<"\n";
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

void do_band_structure_analysis(std::shared_ptr<IOMethods::ResourceHandler> resource)
{
	auto const & opt = resource->get_optns();
	std::shared_ptr<ElectronicBands> interpolSubsetBands;
	// load interpolated bands only if needed.
	auto load_interpol_subset_bands = [&] (
			std::shared_ptr<const ElectronicBands> bands,
			std::vector<int> & bandIndicesInWindow,
			int &firstBandInWindow,
			int &FirstBandNotInWindow) {
		if (interpolSubsetBands)
			return;
		std::vector<double> eneWin = opt.get_ewinbnd();
		if ( eneWin.empty() )
		{
			auto mm = bands->get_min_max();
			eneWin = std::vector<double>{mm.first, mm.second};
		}
		bandIndicesInWindow = bands->get_bands_crossing_energy_window( eneWin );
		auto mm = std::minmax_element(bandIndicesInWindow.begin(), bandIndicesInWindow.end() );
		firstBandInWindow = *(mm.first);
		FirstBandNotInWindow = *(mm.second) + 1;
		auto subsetBands = bands->fft_interpolate_part(
				firstBandInWindow, FirstBandNotInWindow,
				opt.get_fftd(), opt.get_ffts());
		interpolSubsetBands = std::make_shared<ElectronicBands>(std::move(subsetBands));
		};

	int firstBandIn, firstBandNotIn;
	std::vector<int> bandIndicesInEnergyWindow;
	if ( not opt.get_f_mtens().empty() )
	{
		// Find the locations of extrema and print the mass tensor if desired ...
		auto bands = resource->get_electronic_bands_obj();
		load_interpol_subset_bands(bands, bandIndicesInEnergyWindow, firstBandIn, firstBandNotIn);

		int dmeth = 0;
		if ( opt.get_dmeth().compare("fft") == 0 )
			dmeth = 1;
		else if ( opt.get_dmeth().compare("pol") != 0 )
			throw std::runtime_error("unrecognized derivative method");
		write_mass_tensor_file(opt.get_f_mtens(), *interpolSubsetBands, firstBandIn, opt.get_ewinbnd(), dmeth);

		if ( bands->get_bands_crossing_energy_lvls({0.0}).empty() )
		{
			// insulator report vb min and cb max
			output_valenceBandMaxima_conductionBandMinima(boost::filesystem::path(opt.get_f_mtens()));
		}
	}

	if ( not opt.get_f_dos().empty() )
	{
		auto bands = resource->get_dense_electronic_bands_obj();
		load_interpol_subset_bands(bands, bandIndicesInEnergyWindow, firstBandIn, firstBandNotIn);
		auto energySamples = interpolSubsetBands->setup_frequency_grid(
				opt.get_ewinbnd(),
				opt.get_edosnpts());
		interpolSubsetBands->write_tetrahedra_dos_file(opt.get_f_dos(), energySamples);
	}

	if ( not opt.get_f_bands().empty() )
	{
		auto bandsObj = resource->get_dense_electronic_bands_obj();
		auto kpath = resource->get_k_path();

		std::vector<double> bandsAlongPath;
		int numBandsInWindow;
		bandsObj->compute_bands_along_path(
				kpath->get_k_points(),
				resource->get_optns().get_ewinbnd(),
				bandsAlongPath,
				numBandsInWindow,
				resource->get_interpol_reci_tetra_mesh_obj());

		auto gnuplotFile = opt.get_f_bands()+".gp";
		kpath->produce_gnuplot_script_stable_particle(
				gnuplotFile,
				opt.get_f_bands(),
				"\\varepsilon(\\bf{k})",
				bandsAlongPath,
				numBandsInWindow,
				bandsObj->interpret_range(resource->get_optns().get_ewinbnd()));
	}
}

} /* namespace BandStructureAnalysis */
} /* namespace ElectronicStructure */
} /* namespace elephon */
