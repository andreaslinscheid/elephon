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
#include <memory>

namespace elephon
{

// forward declared
namespace LatticeStructure {
	class UnitCell;
	class AtomDisplacement;
	class AtomDisplacementCollection;
	class TetrahedraGrid;
	class PrimitiveToSupercellConnection;
};

namespace ElectronicStructure {
	class ElectronicBands;
	class Wavefunctions;
	class TetrahedraIsosurface;
};

namespace PhononStructure {
	class Phonon;
	class PhononGrid;
	class ForceConstantMatrix;
	class DisplacementPotential;
	class AlphaSquaredF;
};

namespace IOMethods
{

	class KPath;

/**
 * This class is the main interface to load to the required entities in the code only when/if needed.
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

	/**
	 * Construct up the ResourceHandler without loaded data.
	 *
	 * @param dataLoader	The abstracted electronic structure interface that knows how to load specific
	 * 						data dependent to the particular code.
	 */
	ResourceHandler( std::shared_ptr<ElectronicStructureCodeInterface> dataLoader);

	/**
	 * Get access to ElectronicStructureCodeInterface internally stored.
	 * @return	A copy of the pointer to the ElectronicStructureCodeInterface.
	 */
	std::shared_ptr<ElectronicStructureCodeInterface> get_electronic_structure_interface();

	/**
	 * Directly access the input options.
	 *
	 * This is equivalent to ElectronicStructureCodeInterface::get_optns()
	 *
	 * @return	A constant reference to the input options.
	 */
	IOMethods::InputOptions const & get_optns() const;

	/**
	 * obtain the phonon object.
	 *
	 * The object will be initialized to the data as specified via the input options.
	 *
	 * @return Pointer to a constant phonon object.
	 */
	std::shared_ptr<const PhononStructure::Phonon> get_phonon_obj();

	/**
	 * obtain phonon data on a regular grid.
	 *
	 * The object will be initialized to the data as specified via the input options.
	 *
	 * @return Pointer to a constant phonon object.
	 */
	std::shared_ptr<const PhononStructure::PhononGrid> get_phonon_grid_obj();

	std::shared_ptr<const PhononStructure::ForceConstantMatrix> get_forceConstant_obj();

	std::shared_ptr<const LatticeStructure::UnitCell > get_primitive_unitcell_obj();

	std::shared_ptr<const LatticeStructure::UnitCell > get_supercell_obj();

	std::shared_ptr<const LatticeStructure::PrimitiveToSupercellConnection > get_primitive_supercell_connect_obj();

	std::shared_ptr<const LatticeStructure::AtomDisplacementCollection> get_displmts_collection_obj();

	std::shared_ptr<const ElectronicStructure::ElectronicBands> get_electronic_bands_obj();

	/**
	 * obtain the dense electronic band object.
	 *
	 * The code will first try to obtain the 'dense' electronic structure specified in
	 * the 'eld' input option. If that directory does not exist, we gracefully try to load
	 * the band structure in the root directory. Only if these do not exist either, we bail.
	 * Thus the logic is 'get the most dense band structure you can get'.
	 *
	 * @return shared pointer, initialized with the 'dense' band structure.
	 */
	std::shared_ptr<const ElectronicStructure::ElectronicBands> get_dense_electronic_bands_obj();

	std::shared_ptr<const PhononStructure::DisplacementPotential> get_displacement_potential_obj();

	std::shared_ptr<const LatticeStructure::RegularBareGrid> get_real_space_grid_unitcell_obj();

	/**
	 * obtain the wavefunctions object.
	 *
	 * The code will first try to obtain the 'dense' electronic wavefunctions specified in
	 * the 'eld' input option. If that directory does not exist, we gracefully try to load
	 * the wavefunctions in the root directory. Only if these do not exist either, we bail.
	 *
	 * @return shared pointer, initialized with the wavefunction to be used in the code.
	 */
	std::shared_ptr<const ElectronicStructure::Wavefunctions> get_wfct_obj();

	/**
	 * obtain a reciprocal bare mesh according to the input fft dimension.
	 *
	 * @return shared pointer, initialized with a regular bare grid describing the fft mesh.
	 */
	std::shared_ptr<const LatticeStructure::RegularBareGrid> get_interpol_reci_mesh_obj();

	/**
	 * obtain a reciprocal symmetric mesh of tetrahedra according to the input fft dimension.
	 *
	 * @return shared pointer, initialized with a TetrahedraGrid grid describing the fft mesh.
	 */
	std::shared_ptr<const LatticeStructure::TetrahedraGrid> get_interpol_reci_tetra_mesh_obj();

	/**
	 * obtain the tetrahedra mesh of the dense electronic band structure.
	 *
	 * The code will first try to obtain the 'dense' electronic wavefunctions specified in
	 * the 'eld' input option. If that directory does not exist, we gracefully try to load
	 * the wavefunctions in the root directory. Only if these do not exist either, we bail.
	 *
	 * @return shared pointer, initialized with the dense electron band mesh.
	 */
	std::shared_ptr<const LatticeStructure::TetrahedraGrid> get_tetrahedra_grid();

	/**
	 * obtain the tetrahedra iso-surface of the dense electronic band structure.
	 *
	 * Values for the isosurface are taken from the input parameter ea2F.
	 * The code will first try to obtain the 'dense' electronic energies specified in
	 * the 'eld' input option. If that directory does not exist, we gracefully try to load
	 * the electronic structure in the root directory. Only if these do not exist either, we bail.
	 *
	 * @return shared pointer, initialized with the dense electron band mesh.
	 */
	std::shared_ptr<const ElectronicStructure::TetrahedraIsosurface> get_tetrahedra_isosurface();

	/**
	 * obtain the tetrahedra iso-surface of the fft interpolated electronic band structure.
	 *
	 * Values for the isosurface are taken from the input parameter ea2F.
	 * The code will first try to obtain the 'dense' electronic structure specified in
	 * the 'eld' input option. If that directory does not exist, we gracefully try to load
	 * the electronic structure in the root directory. Only if these do not exist either, we bail.
	 * The resulting electronic structure will be interpolated to the fft dim using the fft type interpolation.
	 *
	 * @return shared pointer, initialized with the fft interpolated electron band mesh.
	 */
	std::shared_ptr<const ElectronicStructure::TetrahedraIsosurface> get_tetrahedra_isosurface_fft();

	std::shared_ptr<const IOMethods::KPath> get_k_path();
private:

	std::shared_ptr<ElectronicStructureCodeInterface> dataLoader_;

	std::shared_ptr<PhononStructure::Phonon> ph_;

	std::shared_ptr<PhononStructure::PhononGrid> phGrid_;

	std::shared_ptr<PhononStructure::ForceConstantMatrix> fc_;

	std::shared_ptr< LatticeStructure::UnitCell > primUC_;

	std::shared_ptr< LatticeStructure::UnitCell > scUC_;

	std::shared_ptr< LatticeStructure::PrimitiveToSupercellConnection > primitiveSuperCellConnection_;

	std::shared_ptr<LatticeStructure::AtomDisplacementCollection> displColl_;

	std::shared_ptr<ElectronicStructure::ElectronicBands> bands_;

	std::shared_ptr<ElectronicStructure::ElectronicBands> bandsD_;

	std::shared_ptr<PhononStructure::DisplacementPotential> displPot_;

	std::shared_ptr<LatticeStructure::RegularBareGrid> rsGridUC_;

	std::shared_ptr<ElectronicStructure::Wavefunctions> wfct_;

	std::shared_ptr<LatticeStructure::RegularBareGrid> interpolKGrid_;

	std::shared_ptr<LatticeStructure::TetrahedraGrid> interpolTetraKGrid_;

	std::shared_ptr<LatticeStructure::TetrahedraGrid> tetraGrid_;

	std::shared_ptr<ElectronicStructure::TetrahedraIsosurface> tetraIso_;

	std::shared_ptr<ElectronicStructure::TetrahedraIsosurface> tetraIsoFFT_;

	std::shared_ptr<IOMethods::KPath> kpath_;

	void initialize_phonon_obj();

	void initialize_phonon_grid_obj();

	void initialize_forceConstant_obj();

	void initialize_primitive_unitcell_obj();

	void initialize_supercell_obj();

	void initialize_primitive_supercell_connect_obj();

	void initialize_displmts_collection_obj();

	void initialize_electronic_bands_obj();

	void initialize_dense_electronic_bands_obj();

	void initialize_displacement_potential_obj();

	void initialize_real_space_grid_unitcell_obj();

	void initialize_wfct_obj();

	void initialize_interpol_reci_mesh_obj();

	void initialize_interpol_reci_tetra_mesh_obj();

	void initialize_tetrahedra_grid();

	void initialize_tetrahedra_isosurface();

	void initialize_tetrahedra_isosurface_fft();

	void initialize_k_path();
};

} /* namespace IOMethods */
} /* namespace elephon */

#endif /* ELEPHON_IOMETHODS_RESOURCEHANDLER_H_ */
