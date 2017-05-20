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
#include "IOMethods/InputOptions.h"
#include <vector>
#include <string>
#include <complex>

namespace elephon
{
namespace IOMethods
{

class ElectronicStructureCodeInterface
{
public:

	ElectronicStructureCodeInterface( IOMethods::InputOptions inputOPts );

	virtual ~ElectronicStructureCodeInterface();

	virtual std::vector<std::string> list_all_input_files() const = 0;

	virtual void set_up_run(
			std::string root_directory,
			std::string target_directory,
			std::vector<int> const & kptSampling,
			std::vector<double> const & kptShift,
			LatticeStructure::UnitCell const & unitcell,
			std::map<std::string,std::string> const& options) const = 0;

	virtual std::map<std::string,std::string>
			options_nscf_keep_wfctns_no_relax() const = 0;

	virtual std::map<std::string,std::string>
			options_scf_supercell_no_wfctns_no_relax() const = 0;

	virtual void read_wavefunctions(
			std::vector<std::string> const & files,
			std::vector<int> const & kpts,
			std::vector<int> const & bandIndices,
			std::vector< std::complex<float> > & wfctData,
			std::vector< std::vector<int> > & fourierMap,
			std::vector<int> & fftDim) = 0;

	virtual void read_atoms_list(
			std::vector<std::string> const & baseFiles,
			std::vector<LatticeStructure::Atom> & atoms) = 0;

	virtual void read_cell_paramters(
			std::vector<std::string> const & baseFiles,
			LatticeStructure::LatticeModule & lattice) = 0;

	virtual void read_electronic_potential(
			std::vector<std::string> const & files,
			std::vector<float> & output) = 0;

	virtual void read_symmetries(
			std::vector<std::string> const & files,
			double symPrec,
			LatticeStructure::Symmetry & symmetry) = 0;

	std::vector<std::string> gen_input_file_list( std::string directory ) const;

	virtual void read_kpt_sampling(
			std::string root_directory,
			std::vector<int> & kptSampling,
			std::vector<double> & shifts) = 0;

protected:

	IOMethods::InputOptions inputOPts_;
};

} /* namespace IOMethods */
} /* namespace elephon */

#endif /* ELEPHON_IOMETHODS_ELECTRONICSTRUCTURECODEINTERFACE_H_ */
