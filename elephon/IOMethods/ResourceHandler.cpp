/*	This file ResourceHandler.cpp is part of elephon.
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

#include "IOMethods/ResourceHandler.h"
#include "LatticeStructure/AtomDisplacementCollection.h"
#include "LatticeStructure/PrimitiveToSupercellConnection.h"
#include "IOMethods/KPath.h"
#include "PhononStructure/Phonon.h"
#include "PhononStructure/PhononGrid.h"
#include "PhononStructure/ForceConstantMatrix.h"
#include "LatticeStructure/UnitCell.h"
#include "LatticeStructure/AtomDisplacement.h"
#include "PhononStructure/DisplacementPotential.h"
#include "ElectronicStructure/ElectronicBands.h"
#include "ElectronicStructure/Wavefunctions.h"
#include "ElectronicStructure/TetrahedraIsosurface.h"
#include "PhononStructure/PotentialChangeIrredDisplacement.h"
#include "PhononStructure/Forces.h"
#include <boost/filesystem.hpp>

namespace elephon
{
namespace IOMethods
{

ResourceHandler::ResourceHandler( std::shared_ptr<ElectronicStructureCodeInterface> dataLoader)
	: dataLoader_(std::move(dataLoader))
{

}

std::shared_ptr<ElectronicStructureCodeInterface>
ResourceHandler::get_electronic_structure_interface()
{
	return dataLoader_;
}

IOMethods::InputOptions const &
ResourceHandler::get_optns() const
{
	return dataLoader_->get_optns();
}

std::shared_ptr<const PhononStructure::Phonon>
ResourceHandler::get_phonon_obj()
{
	if ( ! ph_ )
		this->initialize_phonon_obj();
	assert(ph_);
	return ph_;
}

std::shared_ptr<const PhononStructure::PhononGrid>
ResourceHandler::get_phonon_grid_obj()
{
	if ( ! phGrid_ )
		this->initialize_phonon_grid_obj();
	assert(phGrid_);
	return phGrid_;
}

std::shared_ptr<const PhononStructure::ForceConstantMatrix>
ResourceHandler::get_forceConstant_obj()
{
	if ( ! fc_ )
		this->initialize_forceConstant_obj();
	assert(fc_);
	return fc_;
}

std::shared_ptr<const LatticeStructure::UnitCell >
ResourceHandler::get_primitive_unitcell_obj()
{
	if ( ! primUC_ )
		this->initialize_primitive_unitcell_obj();
	assert(primUC_);
	return primUC_;
}

std::shared_ptr<const LatticeStructure::UnitCell >
ResourceHandler::get_supercell_obj()
{
	if ( ! scUC_ )
		this->initialize_supercell_obj();
	assert(scUC_);
	return scUC_;
}

std::shared_ptr<const LatticeStructure::PrimitiveToSupercellConnection >
ResourceHandler::get_primitive_supercell_connect_obj()
{
	if ( ! primitiveSuperCellConnection_ )
		this->initialize_primitive_supercell_connect_obj();
	assert(primitiveSuperCellConnection_);
	return primitiveSuperCellConnection_;
}

std::shared_ptr<const LatticeStructure::AtomDisplacementCollection>
ResourceHandler::get_displmts_collection_obj()
{
	if ( ! displColl_ )
		this->initialize_displmts_collection_obj();
	assert(displColl_);
	return displColl_;
}

std::shared_ptr<const ElectronicStructure::ElectronicBands>
ResourceHandler::get_electronic_bands_obj()
{
	if ( ! bands_ )
		this->initialize_electronic_bands_obj();
	assert(bands_);
	return bands_;
}

std::shared_ptr<const ElectronicStructure::ElectronicBands>
ResourceHandler::get_dense_electronic_bands_obj()
{
	if ( ! bandsD_ )
		this->initialize_dense_electronic_bands_obj();
	assert(bandsD_);
	return bandsD_;
}

std::shared_ptr<const PhononStructure::DisplacementPotential>
ResourceHandler::get_displacement_potential_obj()
{
	if ( ! displPot_ )
		this->initialize_displacement_potential_obj();
	assert(displPot_);
	return displPot_;
}

std::shared_ptr<const LatticeStructure::RegularBareGrid>
ResourceHandler::get_real_space_grid_unitcell_obj()
{
	if ( ! rsGridUC_ )
		this->initialize_real_space_grid_unitcell_obj();
	assert(rsGridUC_);
	return rsGridUC_;
}

std::shared_ptr<const ElectronicStructure::Wavefunctions>
ResourceHandler::get_wfct_obj()
{
	if ( ! wfct_ )
		this->initialize_wfct_obj();
	assert(wfct_);
	return wfct_;
}

std::shared_ptr<const LatticeStructure::RegularBareGrid>
ResourceHandler::get_interpol_reci_mesh_obj()
{
	if ( ! interpolKGrid_ )
		this->initialize_interpol_reci_mesh_obj();
	assert(interpolKGrid_);
	return interpolKGrid_;
}

std::shared_ptr<const LatticeStructure::TetrahedraGrid>
ResourceHandler::get_interpol_reci_tetra_mesh_obj()
{
	if ( ! interpolTetraKGrid_ )
		this->initialize_interpol_reci_tetra_mesh_obj();
	assert(interpolTetraKGrid_);
	return interpolTetraKGrid_;
}

std::shared_ptr<const LatticeStructure::TetrahedraGrid>
ResourceHandler::get_tetrahedra_grid()
{
	if ( ! tetraGrid_ )
		this->initialize_tetrahedra_grid();
	assert(tetraGrid_);
	return tetraGrid_;
}

std::shared_ptr<const ElectronicStructure::TetrahedraIsosurface>
ResourceHandler::get_tetrahedra_isosurface()
{
	if ( ! tetraIso_ )
		this->initialize_tetrahedra_isosurface();
	assert(tetraIso_);
	return tetraIso_;
}

std::shared_ptr<const ElectronicStructure::TetrahedraIsosurface>
ResourceHandler::get_tetrahedra_isosurface_fft()
{
	if ( ! tetraIsoFFT_ )
		this->initialize_tetrahedra_isosurface_fft();
	assert(tetraIsoFFT_);
	return tetraIsoFFT_;
}

std::shared_ptr<const IOMethods::KPath>
ResourceHandler::get_k_path()
{
	if ( ! kpath_ )
		this->initialize_k_path();
	assert(kpath_);
	return kpath_;
}

void
ResourceHandler::initialize_phonon_obj()
{
	auto primUC = this->get_primitive_unitcell_obj();
	std::vector<double> atomicMasses;
	atomicMasses.reserve(primUC->get_atoms_list().size());
	for ( auto it = primUC->get_atoms_list().begin() ; it != primUC->get_atoms_list().end() ; ++it )
		atomicMasses.push_back( it->get_mass() );
	ph_ = std::make_shared<PhononStructure::Phonon>();
	ph_->initialize(this->get_forceConstant_obj(), std::move(atomicMasses) );
}

void
ResourceHandler::initialize_phonon_grid_obj()
{
	PhononStructure::PhononGrid phGrid;
	auto fftgrid = this->get_interpol_reci_mesh_obj();
	auto kgrid = this->get_electronic_bands_obj()->get_grid();
	LatticeStructure::RegularSymmetricGrid symmetricGrid;
	symmetricGrid.initialize( 	fftgrid->get_grid_dim(),
								fftgrid->get_grid_prec(),
								fftgrid->get_grid_shift(),
								kgrid.get_symmetry(),
								kgrid.get_lattice());
	phGrid.initialize(0.0, *this->get_phonon_obj(), std::move(symmetricGrid));
	phGrid_ = std::make_shared<PhononStructure::PhononGrid>(std::move(phGrid));

}

void
ResourceHandler::initialize_forceConstant_obj()
{
	auto phononDir = boost::filesystem::path( dataLoader_->get_optns().get_elphd());

	//Here, we read the forces from the calculator output
	auto forces = std::make_shared<PhononStructure::Forces>();
	forces->initialize(
			this->get_displmts_collection_obj(),
			dataLoader_);

	fc_ = std::make_shared<PhononStructure::ForceConstantMatrix>();
	fc_->initialize(
			this->get_primitive_unitcell_obj(),
			this->get_supercell_obj(),
			this->get_displmts_collection_obj(),
			this->get_primitive_supercell_connect_obj(),
			forces);
}

void
ResourceHandler::initialize_primitive_unitcell_obj()
{
	LatticeStructure::UnitCell uc;
	dataLoader_->read_unit_cell( dataLoader_->get_optns().get_root_dir(), dataLoader_->get_optns().get_gPrec(), uc);
	primUC_ = std::make_shared<LatticeStructure::UnitCell>(std::move(uc));
}

void
ResourceHandler::initialize_supercell_obj()
{
	auto uc = this->get_primitive_unitcell_obj();
	auto scDim = dataLoader_->get_optns().get_scell();
	assert(scDim.size() == 3);
	auto supercell = uc->build_supercell(scDim[0], scDim[1], scDim[2]);
	scUC_ = std::make_shared<LatticeStructure::UnitCell>( std::move(supercell) );
}

void
ResourceHandler::initialize_primitive_supercell_connect_obj()
{
	primitiveSuperCellConnection_ = std::make_shared<LatticeStructure::PrimitiveToSupercellConnection>();
	primitiveSuperCellConnection_->initialize( this->get_primitive_unitcell_obj(), this->get_supercell_obj());
}

void
ResourceHandler::initialize_displmts_collection_obj()
{
	displColl_ = std::make_shared<LatticeStructure::AtomDisplacementCollection>();
	displColl_->initialize( this->get_primitive_unitcell_obj(),
			dataLoader_->get_optns().get_symDispl(),
			dataLoader_->get_optns().get_magdispl());
}

void
ResourceHandler::initialize_electronic_bands_obj()
{
	auto elDir = dataLoader_->get_optns().get_root_dir();
	ElectronicStructure::ElectronicBands bands;
	dataLoader_->read_band_structure(elDir, bands);
	bands_ = std::make_shared<ElectronicStructure::ElectronicBands>(std::move(bands));
}

void
ResourceHandler::initialize_dense_electronic_bands_obj()
{
	auto elDir = dataLoader_->get_optns().get_eld();

	// if the path 'eld' does not exist, we gracefully try
	// the location of the root dir. Only neither one works, we bail.
	if ( not elDir.empty() )
	{
		boost::filesystem::path eldense_path(elDir);
		if ( not boost::filesystem::exists(eldense_path) )
			elDir.clear();
	}

	if ( elDir.empty() )
	{
		auto normalBands = this->get_electronic_bands_obj();
		bandsD_ = std::make_shared<ElectronicStructure::ElectronicBands>(std::move(*normalBands));
	}
	else
	{
		ElectronicStructure::ElectronicBands bands;
		dataLoader_->read_band_structure(elDir, bands);
		bandsD_ = std::make_shared<ElectronicStructure::ElectronicBands>(std::move(bands));
	}
}

void
ResourceHandler::initialize_displacement_potential_obj()
{
	auto supercell = this->get_supercell_obj();
	auto primitiveCell = this->get_primitive_unitcell_obj();
	auto phononDir = boost::filesystem::path( dataLoader_->get_optns().get_elphd());

	const double gPrec = dataLoader_->get_optns().get_gPrec();
	auto displColl = this->get_displmts_collection_obj();
	const int nIrdDispl = displColl->get_tota_num_irred_displacements();
	if ( nIrdDispl < 1)
		throw std::runtime_error("Problem getting number of irreducible displacements");

	// here we load the grid data for both supercell and primitive cell and their grids
	auto rsGridSC = std::make_shared<LatticeStructure::RegularBareGrid>();

	// Read the normal periodic potential form the root dir
	LatticeStructure::DataRegularAndRadialGrid<double> primitiveCellPotential;
	std::vector<int> dimPCGrid;
	dataLoader_->read_electronic_potential(
			dataLoader_->get_optns().get_root_dir(),
			dimPCGrid,
			primitiveCellPotential);
	auto rsGridPC = std::make_shared<LatticeStructure::RegularBareGrid>();
	rsGridPC->initialize( dimPCGrid, false, gPrec, {0.0, 0.0, 0.0}, primitiveCell->get_lattice());

	// Here, we read the potential from the calculator output
	std::vector<std::shared_ptr<const PhononStructure::PotentialChangeIrredDisplacement>> displPot( nIrdDispl );
	int idispl = 0;
	for ( auto atomDispl : displColl->get_irreducible_displacements() )
	{
		for ( auto const & irreducibleDispl : atomDispl.second )
		{
			std::vector<int> dimSCGrid;

			// regular grid part
			LatticeStructure::DataRegularAndRadialGrid<double> thisDisplPot;
			dataLoader_->read_electronic_potential(
					( phononDir / ("displ_"+std::to_string(idispl)) ).string(),
					dimSCGrid,
					thisDisplPot);

			if (idispl == 0 )
				rsGridSC->initialize( dimSCGrid, false, gPrec, {0.0, 0.0, 0.0}, supercell->get_lattice());

			for (int ix = 0 ; ix < 3 ; ++ix)
				if ( dimSCGrid[ix] % rsGridSC->get_grid_dim()[ix] )
					throw std::runtime_error("Input error: The Fourier grid of the supercell is not consistent among displacements.");


			PhononStructure::PotentialChangeIrredDisplacement potChangeDispl;
			potChangeDispl.initialize(	irreducibleDispl,
										primitiveCellPotential,
										thisDisplPot,
										rsGridPC,
										rsGridSC,
										this->get_primitive_supercell_connect_obj());

			displPot[idispl] = std::make_shared<PhononStructure::PotentialChangeIrredDisplacement>(std::move(potChangeDispl));
			++idispl;
		}
	}

	// check if the grids are compatible
	for (int ix = 0 ; ix < 3 ; ++ix)
		if ( rsGridSC->get_grid_dim()[ix] % rsGridPC->get_grid_dim()[ix] )
			throw std::runtime_error("Input error: The Fourier grid of the supercell is not a multiple of the unit cell grid.");

	displPot_ = std::make_shared<PhononStructure::DisplacementPotential>();
	displPot_->initialize(
			this->get_primitive_unitcell_obj(),
			this->get_supercell_obj(),
			displColl,
			this->get_primitive_supercell_connect_obj(),
			*rsGridPC,
			*rsGridSC,
			displPot);
}

void
ResourceHandler::initialize_real_space_grid_unitcell_obj()
{
	auto root_dir = boost::filesystem::path( dataLoader_->get_optns().get_root_dir());
	auto gPrec = dataLoader_->get_optns().get_gPrec();

	auto unitcell = this->get_primitive_unitcell_obj();

	rsGridUC_ = std::make_shared<LatticeStructure::RegularBareGrid>();
	rsGridUC_->initialize(dataLoader_->read_charge_real_space_grid_dim(root_dir.string()),
						  /* is reciprocal mesh */false,
						  gPrec,
						  {0.0, 0.0, 0.0},
						  unitcell->get_lattice() );
}

void
ResourceHandler::initialize_wfct_obj()
{
	ElectronicStructure::Wavefunctions wfcts;
	std::string wfctdir = dataLoader_->get_optns().get_eld();

	// if the path 'eld' does not exist, we gracefully try
	// the location of the root dir. Only neither one works, we bail.
	if ( not wfctdir.empty() )
	{
		boost::filesystem::path wfct_path(wfctdir);
		if ( not boost::filesystem::exists(wfct_path) )
			wfctdir.clear();
	}

	if ( wfctdir.empty() )
		wfctdir = dataLoader_->get_optns().get_root_dir();

	auto denseBands = this->get_dense_electronic_bands_obj();

	wfcts.initialize(wfctdir, denseBands->get_grid(), denseBands->get_nBnd(), dataLoader_ );
	wfct_ = std::make_shared<ElectronicStructure::Wavefunctions>(std::move(wfcts));
}

void
ResourceHandler::initialize_interpol_reci_mesh_obj()
{
	// the not-interpolated k grid is contained in the band structure
	auto bands = this->get_electronic_bands_obj();
	auto interpolDim = bands->get_grid().interpret_fft_dim_input( dataLoader_->get_optns().get_fftd() );

	interpolKGrid_ = std::make_shared<LatticeStructure::RegularBareGrid>();
	interpolKGrid_->initialize(
			interpolDim,
			true,
			bands->get_grid().get_grid_prec(),
			dataLoader_->get_optns().get_ffts(),
			bands->get_grid().get_lattice());
}

void
ResourceHandler::initialize_interpol_reci_tetra_mesh_obj()
{
	auto symmetry = this->get_primitive_unitcell_obj()->get_symmetry();
	symmetry.set_reciprocal_space_sym(true);
	auto bareGrid = this->get_interpol_reci_mesh_obj();
	LatticeStructure::RegularSymmetricGrid symmetricGrid;
	symmetricGrid.initialize( 	bareGrid->get_grid_dim(),
								bareGrid->get_grid_prec(),
								bareGrid->get_grid_shift(),
								symmetry,
								bareGrid->get_lattice());

	interpolTetraKGrid_ = std::make_shared<LatticeStructure::TetrahedraGrid>();
	interpolTetraKGrid_->initialize(std::make_shared<LatticeStructure::RegularSymmetricGrid>(std::move(symmetricGrid)));
}

void
ResourceHandler::initialize_tetrahedra_grid()
{
	LatticeStructure::TetrahedraGrid tg;
	tg.initialize( std::make_shared<elephon::LatticeStructure::RegularSymmetricGrid>(
			this->get_dense_electronic_bands_obj()->get_grid() ));
	tetraGrid_ = std::make_shared<LatticeStructure::TetrahedraGrid>( std::move(tg) );
}

void
ResourceHandler::initialize_tetrahedra_isosurface()
{
	ElectronicStructure::TetrahedraIsosurface tis;
	tis.initialize( this->get_tetrahedra_grid(),
					this->get_dense_electronic_bands_obj(),
					this->get_optns().get_ea2f() );
	tetraIso_ = std::make_shared<ElectronicStructure::TetrahedraIsosurface>( std::move(tis) );
}

void
ResourceHandler::initialize_tetrahedra_isosurface_fft()
{
	auto fft_tetra_mesh = this->get_interpol_reci_tetra_mesh_obj();
	ElectronicStructure::ElectronicBands denseBands = *this->get_dense_electronic_bands_obj();
	auto bareGrid = this->get_interpol_reci_mesh_obj();
	denseBands.fft_interpolate(bareGrid->get_grid_dim(), bareGrid->get_grid_shift());
	auto fftbands_ptr = std::make_shared<ElectronicStructure::ElectronicBands>(std::move(denseBands));
	ElectronicStructure::TetrahedraIsosurface tis;
	tis.initialize( this->get_interpol_reci_tetra_mesh_obj(),
					fftbands_ptr,
					this->get_optns().get_ea2f() );
	tetraIsoFFT_ = std::make_shared<ElectronicStructure::TetrahedraIsosurface>( std::move(tis) );
}

void
ResourceHandler::initialize_k_path()
{
	kpath_ = std::make_shared<IOMethods::KPath>();
	kpath_->read_kpath_file(this->get_optns().get_f_kpath());
}

} /* namespace IOMethods */
} /* namespace elephon */
