/*	This file ResourceHandler.h is part of elephon.
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
 *  Created on: Oct 4, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_IOMETHODS_RESOURCEHANDLER_H_
#define ELEPHON_IOMETHODS_RESOURCEHANDLER_H_

#include "IOMethods/ElectronicStructureCodeInterface.h"
#include "PhononStructure/Phonon.h"
#include "PhononStructure/ForceConstantMatrix.h"
#include "LatticeStructure/UnitCell.h"
#include "LatticeStructure/AtomDisplacement.h"
#include "PhononStructure/DisplacementPotential.h"
#include "ElectronicStructure/ElectronicBands.h"
#include "ElectronicStructure/Wavefunctions.h"
#include "ElectronicStructure/TetrahedraIsosurface.h"
#include <memory>

namespace elephon
{
namespace IOMethods
{

/**
 * This class is the main interface to the required entities in the code.
 *
 * It stores the phonon, force constant, displacement potential, electronic grid and so on
 * if they have been computed and triggers their calculation if they are needed. The logic is
 * that only one variant of each of these resources is needed. If they locally need to be overwritten
 * this will have to be treated locally. Thus, this method does not provide write access to resources
 * pointed to.
 */
class ResourceHandler
{
public:

	ResourceHandler( std::shared_ptr<ElectronicStructureCodeInterface> dataLoader);

	IOMethods::InputOptions const & get_optns() const;

	std::shared_ptr<const PhononStructure::Phonon> get_phonon_obj();

	std::shared_ptr<const PhononStructure::ForceConstantMatrix> get_forceConstant_obj();

	std::shared_ptr<const LatticeStructure::UnitCell > get_primitive_unitcell_obj();

	std::shared_ptr<const LatticeStructure::UnitCell > get_supercell_obj();

	std::shared_ptr<const std::vector<LatticeStructure::AtomDisplacement> > get_irrd_displmts_obj();

	std::shared_ptr<const ElectronicStructure::ElectronicBands> get_electronic_bands_obj();

	std::shared_ptr<const ElectronicStructure::ElectronicBands> get_dense_electronic_bands_obj();

	std::shared_ptr<const PhononStructure::DisplacementPotential> get_displacement_potential_obj();

	std::shared_ptr<const ElectronicStructure::Wavefunctions> get_wfct_obj();

	std::shared_ptr<const LatticeStructure::RegularBareGrid> get_interpol_reci_mesh_obj();

	/**
	 * obtain the tetrahedra mesh of the dense electronic band structure.
	 *
	 * @return shared pointer, initialized with the dense electron band mesh.
	 */
	std::shared_ptr<const LatticeStructure::TetrahedraGrid> get_tetrahedra_grid();

	/**
	 * obtain the tetrahedra iso-surface of the dense electronic band structure.
	 *
	 * Values for the isosurface are taken from the input parameter ea2F.
	 *
	 * @return shared pointer, initialized with the dense electron band mesh.
	 */
	std::shared_ptr<const ElectronicStructure::TetrahedraIsosurface> get_tetrahedra_isosurface();

private:

	std::shared_ptr<ElectronicStructureCodeInterface> dataLoader_;

	std::shared_ptr<PhononStructure::Phonon> ph_;

	std::shared_ptr<PhononStructure::ForceConstantMatrix> fc_;

	std::shared_ptr< LatticeStructure::UnitCell > primUC_;

	std::shared_ptr< LatticeStructure::UnitCell > scUC_;

	std::shared_ptr<std::vector<LatticeStructure::AtomDisplacement> > irredDispl_;

	std::shared_ptr<ElectronicStructure::ElectronicBands> bands_;

	std::shared_ptr<ElectronicStructure::ElectronicBands> bandsD_;

	std::shared_ptr<PhononStructure::DisplacementPotential> displPot_;

	std::shared_ptr<ElectronicStructure::Wavefunctions> wfct_;

	std::shared_ptr<LatticeStructure::RegularBareGrid> interpolKGrid_;

	std::shared_ptr<LatticeStructure::TetrahedraGrid> tetraGrid_;

	std::shared_ptr<ElectronicStructure::TetrahedraIsosurface> tetraIso_;

	void initialize_phonon_obj();

	void initialize_forceConstant_obj();

	void initialize_primitive_unitcell_obj();

	void initialize_supercell_obj();

	void initialize_irrd_displmts_obj();

	void initialize_electronic_bands_obj();

	void initialize_dense_electronic_bands_obj();

	void initialize_displacement_potential_obj();

	void initialize_wfct_obj();

	void initialize_interpol_reci_mesh_obj();

	void initialize_tetrahedra_grid();

	void initialize_tetrahedra_isosurface();
};

} /* namespace IOMethods */
} /* namespace elephon */

#endif /* ELEPHON_IOMETHODS_RESOURCEHANDLER_H_ */
