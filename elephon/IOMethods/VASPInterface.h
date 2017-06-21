/*	This file VASPInterface.h is part of elephon.
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
 *  Created on: May 17, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_IOMETHODS_VASPINTERFACE_H_
#define ELEPHON_IOMETHODS_VASPINTERFACE_H_

#include "IOMethods/ElectronicStructureCodeInterface.h"
#include "IOMethods/ReadVASPPoscar.h"
#include "IOMethods/ReadVASPSymmetries.h"
#include "IOMethods/ReadVASPWaveFunction.h"
#include "IOMethods/ReadVASPPotcar.h"
#include "IOMethods/ReadVASPxmlFile.h"
#include <vector>
#include <map>

namespace elephon
{
namespace IOMethods
{

class VASPInterface : public ElectronicStructureCodeInterface
{
public:

	using ElectronicStructureCodeInterface::ElectronicStructureCodeInterface;

	void set_up_run(
			std::string root_directory,
			std::string target_directory,
			std::vector<int> const & kptSampling,
			std::vector<double> const & kptShift,
			LatticeStructure::UnitCell const & unitcell,
			std::map<std::string,std::string> const& options) const;

	std::map<std::string,std::string>
			options_nscf_keep_wfctns_no_relax() const;

	std::map<std::string,std::string>
			options_scf_supercell_no_wfctns_no_relax() const;

	void read_wavefunctions(
			std::string root_directory,
			std::vector<int> const & kpts,
			std::vector<int> const & bandIndices,
			std::vector< std::complex<float> > & wfctData,
			std::vector< int > & npwPerKpt);

	std::vector<int> get_max_fft_dims() const;

	void compute_fourier_map(
			std::vector<double> const & kpts,
			std::vector< std::vector<int> > & fourierMap);

	void read_cell_paramters(
			std::string root_directory,
			double symPrec,
			LatticeStructure::RegularGrid & kPointMesh,
			LatticeStructure::LatticeModule & lattice,
			std::vector<LatticeStructure::Atom> & atoms,
			LatticeStructure::Symmetry & symmetry);

	void read_forces(
			std::string root_directory,
			std::vector<double> & forces);

	void read_electronic_potential(
			std::string root_directory,
			std::vector<int> & dims,
			std::vector<double> & output);

	void read_kpt_sampling(
			std::string root_directory,
			std::vector<int> & kptSampling,
			std::vector<double> & shifts) const;

	void read_electronic_structure(
			std::string root_directory,
			int & nBnd,
			int & nkptsIrred,
			std::vector<double> & energies,
			double & fermiEnergy);

private:

	ReadVASPPoscar posReader_;

	ReadVASPPotcar potReader_;

	ReadVASPSymmetries symReader_;

	ReadVASPWaveFunction wfcReader_;

	ReadVASPxmlFile xmlReader_;

	void overwrite_POSCAR_file( std::string filename,
			std::vector<std::string > const & potcarAtomOrder,
			LatticeStructure::UnitCell const & unitcell ) const ;

	void write_KPOINTS_file(std::string filename,
			std::vector<double> const & kptShift,
			std::vector<int> const & monkhPackGrid) const;

	void modify_incar_file(std::string filename,
			std::map<std::string,std::string> const & optionsToBeReset) const;

	std::vector<std::string > read_potcar_atom_order( std::string filename ) const;

	std::string get_textfile_constent( std::string filename ) const;

	void check_open_poscar(std::string const & root_dir );
};

} /* namespace IOMethods */
} /* namespace elephon */

#endif /* ELEPHON_IOMETHODS_VASPINTERFACE_H_ */
