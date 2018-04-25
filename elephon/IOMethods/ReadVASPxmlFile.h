/*	This file ReadVASPxmlFile.h is part of elephon.
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
 *  Created on: May 31, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_IOMETHODS_READVASPXMLFILE_H_
#define ELEPHON_IOMETHODS_READVASPXMLFILE_H_

#include "LatticeStructure/Atom.h"
#include <boost/property_tree/ptree.hpp>
#include <vector>
#include <string>

namespace elephon
{
namespace IOMethods
{

/**
 * This class handles the reading of data from the vasprun.xml file.
 */
class ReadVASPxmlFile
{
public:

	bool is_parsed() const;

	void parse_file( std::string filename );

	std::vector<double> const & get_forces();

	std::vector<double> const & get_k_points();

	double get_Fermi_energy();

	std::vector<double> const & get_energies();

	int get_nBnd();

	int get_nkp();

	std::vector<double> const & get_lattice_matrix();

	std::vector<LatticeStructure::Atom> const & get_atoms_list();

	std::vector<int> get_wfct_fourier_dim();

	std::vector<int> get_charge_fourier_dim();

	std::vector<int> get_k_grid_dim();

	std::vector<double> get_k_grid_shift();
private:

	int nBnd_ = 0;

	double eFermi_ = 0;

	std::string filename_;

	bool parseForces_ = true;
	std::vector<double> forces_;

	bool parseKPoints_ = true;
	std::vector<double> kpoints_;

	bool parseEnergies_ = true;
	std::vector<double> energies_;

	bool parseLatticeMatrix_ = true;
	std::vector<double> latticeMat_;

	bool parseAtoms_ = true;
	std::vector<LatticeStructure::Atom> atoms_;

	std::vector<int> kDim_;

	std::vector<double> kShift_;

	bool parseFFTDims_ = true;
	std::vector<int> wfctFourierDim_;

	std::vector<int> chargeFourierDim_;

	boost::property_tree::ptree pt_;

	void parse_forces();

	void parse_latticeMat();

	void parse_atoms();

	void parse_energies();

	void parse_kpoints();

	void parse_fftdims();
};

} /* namespace IOMethods */
} /* namespace elephon */

#endif /* ELEPHON_IOMETHODS_READVASPXMLFILE_H_ */
