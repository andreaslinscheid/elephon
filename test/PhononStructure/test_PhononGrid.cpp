/*	This file test_PhononGrid.cpp is part of elephon.
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
 *  Created on: Oct 3, 2017
 *      Author: A. Linscheid
 */

#define BOOST_TEST_MODULE PhononStructure
#include <boost/test/unit_test.hpp>
#include "fixtures/MockStartup.h"
#include "fixtures/DataLoader.h"
#include "PhononStructure/PhononGrid.h"
#include "fixtures/scenarios.h"

BOOST_AUTO_TEST_CASE( phononDOS_write_file_vasp_fcc_primitive )
{
	auto resourceHandler = elephon::test::fixtures::scenarios::load_Al_fcc_primitive_vasp_sc2x2x2();

	auto bands = resourceHandler->get_electronic_bands_obj();

	elephon::LatticeStructure::RegularSymmetricGrid qGrid;
	qGrid.initialize(
			{20, 20, 20},
			resourceHandler->get_optns().get_gPrec(),
			{0.0, 0.0, 0.0},
			bands->get_grid().get_symmetry(),
			bands->get_grid().get_lattice());

	auto ph = resourceHandler->get_phonon_obj();

	elephon::PhononStructure::PhononGrid phgrid;
	phgrid.initialize(0.0, *ph, qGrid);

	double freqmin = 0;
	double freqmax = 10;
	int numFreq = 100;
	std::vector<double> frequencies(numFreq);
	for ( int iw = 0 ; iw < numFreq; ++iw)
		frequencies[iw] = freqmin + (freqmax - freqmin) * double(iw + 0.5) / (numFreq);

	elephon::LatticeStructure::TetrahedraGrid tetra;
	tetra.initialize( std::make_shared<decltype(qGrid)>(qGrid) );

	auto phDosFilename = boost::filesystem::path(resourceHandler->get_optns().get_f_phdos());
	boost::filesystem::remove(phDosFilename);
	phgrid.write_phonon_dos_file(phDosFilename.string(), frequencies, std::make_shared<decltype(tetra)>(tetra));

	std::ifstream file( phDosFilename.c_str() );
	if ( ! file.good() )
		throw std::runtime_error("Problem opening file");

	std::vector<double> phDos, freq;
	std::string line;
	std::getline(file, line); // comment
	while ( std::getline(file, line) )
	{
		std::size_t sz;
 		freq.push_back( std::stof(line, &sz) );
 		phDos.push_back( std::stof(line.substr(sz)) );
	}
	BOOST_REQUIRE_EQUAL(phDos.size(), numFreq);
	BOOST_CHECK_SMALL( (freq[1]-freq[0]) - (freqmax - freqmin) / double(numFreq) , 1e-6);

	double nM = 0;
	for ( auto phD : phDos )
		nM += phD*(freqmax - freqmin) / double(numFreq);

	std::cout << "Integrated phonon DOS : "<< nM << " with a number of modes of " << ph->get_num_modes() << std::endl;
	BOOST_CHECK_CLOSE( nM , phgrid.get_nData_gpt() , 10);
	boost::filesystem::remove(phDosFilename);

	std::vector<double> phDos_wan;
	phgrid.compute_DOS_wan(*ph, frequencies, phDos_wan);

	std::vector<double> phDos_poly;
	phgrid.compute_DOS(frequencies, phDos_poly);

	BOOST_CHECK_EQUAL( phDos_poly.size() , phDos_wan.size());
	BOOST_CHECK_EQUAL( phDos_poly.size() , phDos.size());
	for (int iw = 0 ; iw < phDos.size() ; ++iw )
		std::cout << frequencies[iw] << '\t' << phDos[iw] << '\t' << phDos_wan[iw] << '\t' << phDos_poly[iw] << '\n';
}



