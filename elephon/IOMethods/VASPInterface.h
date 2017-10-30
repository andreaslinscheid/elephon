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
#include "ElectronicStructure/ElectronicBands.h"
#include "IOMethods/ReadVASPPoscar.h"
#include "IOMethods/ReadVASPSymmetries.h"
#include "IOMethods/ReadVASPWaveFunction.h"
#include "IOMethods/ReadVASPxmlFile.h"
#include "IOMethods/ReadVASPKpoints.h"
#include "LatticeStructure/Atom.h"
#include <vector>
#include <map>
#include "ReadVASPLocpot.h"

namespace elephon
{
namespace IOMethods
{

class VASPInterface : public ElectronicStructureCodeInterface
{
public:

	using ElectronicStructureCodeInterface::ElectronicStructureCodeInterface;

	std::string code_tag() const;

	void check_prep_run(
			std::string root_directory ) const;

	void set_up_run(
			std::string root_directory,
			std::string target_directory,
			std::vector<int> const & kptSampling,
			std::vector<double> const & kptShift,
			LatticeStructure::UnitCell const & unitcell,
			std::map<std::string,std::string> const& options) const;

	void copy_charge(
			std::string root_directory,
			std::string target_directory) const;

	std::map<std::string,std::string>
			options_nscf_keep_wfctns_no_relax() const;

	std::map<std::string,std::string>
			options_scf_supercell_no_wfctns_no_relax() const;

	void read_wavefunctions(
			std::string root_directory,
			std::vector<int> const & kpts,
			std::vector<int> const & bandIndices,
			std::vector< std::vector< std::complex<float> > > & wfctData,
			std::vector< int > & npwPerKpt);

	std::vector<int> get_max_fft_dims();

	void compute_fourier_map(
			std::vector<double> const & kpts,
			std::vector< std::vector<int> > & fourierMap,
			double gridPrec);

	void read_cell_paramters(
			std::string root_directory,
			double symPrec,
			LatticeStructure::RegularSymmetricGrid & kPointMesh,
			LatticeStructure::LatticeModule & lattice,
			std::vector<LatticeStructure::Atom> & atoms,
			LatticeStructure::Symmetry & symmetry);

	std::vector<int>
			read_wfct_real_space_grid_dim(std::string root_directory);

	std::vector<int>
			read_charge_real_space_grid_dim(std::string root_directory);

	void read_unit_cell(
			std::string root_directory,
			double symprec,
			LatticeStructure::UnitCell & unitcell );

	void read_lattice_structure(
			std::string root_directory,
			LatticeStructure::LatticeModule & lattice);

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
			std::vector<double> & shifts);

	void read_reciprocal_symmetric_grid(
			std::string root_directory,
			LatticeStructure::RegularSymmetricGrid & kgrid);

	void read_band_structure(
			std::string root_directory,
			ElectronicStructure::ElectronicBands & bands);

	void read_electronic_structure(
			std::string root_directory,
			int & nBnd,
			int & nkptsIrred,
			std::vector<double> & energies,
			double & fermiEnergy);

	void read_atoms_list(
			std::string root_directory,
			std::vector<LatticeStructure::Atom> & atoms);

	void read_nBnd(
			std::string root_directory,
			int & nBnd);

	/**
	 * read the symmetry of the system from the OUTCAR file.
	 *
	 * NOTE: the symmetry will be initially set to a real-space
	 * symmetry, i.e. time reversal symmetry is not automatically an inversion.
	 *
	 * @param root_directory	The directory where to look for the OUTCAR file.
	 * @param symprec			The precision of the lattice.
	 * @param lattice			The description of the cell.
	 * @param symmetry			the object to be set to the symmetry of the system.
	 */
	void read_symmetry(
			std::string root_directory,
			double symprec,
			LatticeStructure::LatticeModule const& lattice,
			LatticeStructure::Symmetry & symmetry);

private:

	ReadVASPPoscar posReader_;

	ReadVASPLocpot potReader_;

	ReadVASPSymmetries symReader_;

	ReadVASPWaveFunction wfcReader_;

	ReadVASPxmlFile xmlReader_;

	ReadVASPKpoints kpointReader_;

	void overwrite_POSCAR_file( std::string filename,
			std::vector<std::string > const & potcarAtomOrder,
			LatticeStructure::UnitCell const & unitcell ) const ;

	void write_KPOINTS_file(std::string filename,
			std::vector<double> const & kptShift,
			std::vector<int> const & monkhPackGrid) const;

	void modify_incar_file(std::string filename,
			std::map<std::string,std::string> const & optionsToBeReset) const;

	std::vector<std::pair<std::string, double> > read_potcar_atom_order( std::string filename ) const;

	std::string get_textfile_content( std::string filename ) const;

	void check_open_poscar(std::string const & root_dir );
};

} /* namespace IOMethods */
} /* namespace elephon */

#endif /* ELEPHON_IOMETHODS_VASPINTERFACE_H_ */
