/*	This file ElectronicStructureCodeInterface.h is part of elephon.
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

#ifndef ELEPHON_IOMETHODS_ELECTRONICSTRUCTURECODEINTERFACE_H_
#define ELEPHON_IOMETHODS_ELECTRONICSTRUCTURECODEINTERFACE_H_

#include "LatticeStructure/UnitCell.h"
#include "ElectronicStructure/ElectronicBands.h"
#include "IOMethods/InputOptions.h"
#include "LatticeStructure/RegularSymmetricGrid.h"
#include <vector>
#include <string>
#include <complex>
#include <memory>

namespace elephon
{
namespace IOMethods
{

class ElectronicStructureCodeInterface
{
public:

	ElectronicStructureCodeInterface( IOMethods::InputOptions inputOPts );

	IOMethods::InputOptions const & get_optns() const;

	virtual std::string code_tag() const = 0;

	virtual ~ElectronicStructureCodeInterface();

	virtual void check_prep_run(
			std::string root_directory ) const = 0;

	virtual void set_up_run(
			std::string root_directory,
			std::string target_directory,
			std::vector<int> const & kptSampling,
			std::vector<double> const & kptShift,
			LatticeStructure::UnitCell const & unitcell,
			std::map<std::string,std::string> const& options) const = 0;

	virtual void copy_charge(
			std::string root_directory,
			std::string target_directory)  const = 0;

	virtual std::map<std::string,std::string>
			options_nscf_keep_wfctns_no_relax() const = 0;

	virtual std::map<std::string,std::string>
			options_scf_supercell_no_wfctns_no_relax() const = 0;

	virtual std::vector<int> get_max_fft_dims() = 0;

	/**
	 * Obtain the wave function data.
	 *
	 * Sign convention: The code uses the following sign and normalization convention
	 * w(G) =      sum(r) w(G) exp(+iG.r)
	 * w(r) = 1/NG sum(G) w(G) exp(-iG.r)
	 *
	 * @param root_directory
	 * @param kpts
	 * @param bandIndices
	 * @param wfctData
	 * @param npwPerKpt
	 */
	virtual void read_wavefunctions(
			std::string root_directory,
			std::vector<int> const & kpts,
			std::vector<int> const & bandIndices,
			std::vector< std::vector< std::complex<float> > > & wfctData,
			std::vector< int > & npwPerKpt) = 0;

	virtual void compute_fourier_map(
			std::vector<double> const & kpts,
			std::vector< std::vector<int> > & fourierMap,
			double gridPrec) = 0;

	virtual void read_cell_paramters(
			std::string root_directory,
			double symPrec,
			LatticeStructure::RegularSymmetricGrid & kPointMesh,
			LatticeStructure::LatticeModule & lattice,
			std::vector<LatticeStructure::Atom> & atoms,
			LatticeStructure::Symmetry & symmetry) = 0;

	virtual std::vector<int>
			read_wfct_real_space_grid_dim(std::string root_directory) = 0;

	virtual std::vector<int>
			read_charge_real_space_grid_dim(std::string root_directory) = 0;

	virtual void read_unit_cell(
			std::string root_directory,
			double symprec,
			LatticeStructure::UnitCell & unitcell ) = 0;

	virtual void read_lattice_structure(
			std::string root_directory,
			LatticeStructure::LatticeModule & lattice) = 0;

	virtual void read_electronic_potential(
			std::string root_directory,
			std::vector<int> & dims,
			std::vector<double> & output) = 0;

	virtual void read_band_structure(
			std::string root_directory,
			ElectronicStructure::ElectronicBands & bands) = 0;

	virtual void read_reciprocal_symmetric_grid(
			std::string root_directory,
			LatticeStructure::RegularSymmetricGrid & kgrid) = 0;

	//The returned forces must be in the same order as the atoms in the unperturbed supercell!
	virtual void read_forces(
			std::string root_directory,
			std::vector<double> & forces) = 0;

	virtual void read_kpt_sampling(
			std::string root_directory,
			std::vector<int> & kptSampling,
			std::vector<double> & shifts) = 0;

	virtual void read_nBnd(
			std::string root_directory,
			int & nBnd) = 0;

	std::vector<std::string> gen_input_file_list( std::string directory ) const;

private:

	IOMethods::InputOptions inputOPts_;
};

} /* namespace IOMethods */
} /* namespace elephon */

#endif /* ELEPHON_IOMETHODS_ELECTRONICSTRUCTURECODEINTERFACE_H_ */
