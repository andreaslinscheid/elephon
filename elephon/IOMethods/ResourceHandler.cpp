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
#include <boost/filesystem.hpp>

namespace elephon
{
namespace IOMethods
{

ResourceHandler::ResourceHandler( std::shared_ptr<ElectronicStructureCodeInterface> dataLoader)
	: dataLoader_(std::move(dataLoader))
{

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

std::shared_ptr<const std::vector<LatticeStructure::AtomDisplacement> >
ResourceHandler::get_irrd_displmts_obj()
{
	if ( ! irredDispl_ )
		this->initialize_irrd_displmts_obj();
	assert(irredDispl_);
	return irredDispl_;
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
ResourceHandler::initialize_forceConstant_obj()
{
	auto phononDir = boost::filesystem::path( dataLoader_->get_optns().get_elphd());

	//Here, we read the forces from the calculator output
	int nIrdDispl = int( this->get_irrd_displmts_obj()->size());
	std::vector<std::vector<double>> forces( nIrdDispl );
	std::vector<double> thisForces;
	for ( int idispl = 0 ; idispl < nIrdDispl; ++idispl )
	{
		dataLoader_->read_forces(
				(phononDir / (std::string("displ_")+std::to_string(idispl))).string(),
				thisForces);
		forces[idispl] = std::move(thisForces);
	}

	fc_ = std::make_shared<PhononStructure::ForceConstantMatrix>();
	fc_->build(
			this->get_primitive_unitcell_obj(),
			this->get_supercell_obj(),
			this->get_irrd_displmts_obj(),
			std::move(forces));
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
ResourceHandler::initialize_irrd_displmts_obj()
{
	auto ops = dataLoader_->get_optns();
	auto uc = this->get_primitive_unitcell_obj();
	std::vector<LatticeStructure::AtomDisplacement> irredDispl;
	uc->generate_displacements(
			ops.get_magdispl(),
			ops.get_symDispl(),
			irredDispl);
	irredDispl_ = std::make_shared<std::vector<LatticeStructure::AtomDisplacement>>(std::move(irredDispl));
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
	if ( elDir.empty() )
		bandsD_ = bands_;
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
	auto unitCell = this->get_primitive_unitcell_obj();
	auto phononDir = boost::filesystem::path( dataLoader_->get_optns().get_elphd());

	auto gPrec = dataLoader_->get_optns().get_gPrec();

	// Here, we read the potential from the calculator output
	int nIrdDispl = int( this->get_irrd_displmts_obj()->size());
	std::vector<std::vector<double>> displPot( nIrdDispl );
	std::vector<int> dim;
	for ( int idispl = 0 ; idispl < nIrdDispl; ++idispl )
	{
		std::vector<double> thisDisplPot;
		dataLoader_->read_electronic_potential(
				( phononDir / ("displ_"+std::to_string(idispl)) ).string(),
				dim,
				thisDisplPot);

		displPot[idispl] = std::move(thisDisplPot);
	}
	LatticeStructure::RegularBareGrid rsGridSC;
	rsGridSC.initialize( dim, false, gPrec, {0.0, 0.0, 0.0}, supercell->get_lattice());

	//Read the normal periodic potential form the root dir
	std::vector<double> primitiveCellPotential;
	dataLoader_->read_electronic_potential(
			dataLoader_->get_optns().get_root_dir(),
			dim,
			primitiveCellPotential);
	elephon::LatticeStructure::RegularBareGrid rsGridUC;
	rsGridUC.initialize( dim, false, gPrec, {0.0, 0.0, 0.0}, unitCell->get_lattice());

	displPot_ = std::make_shared<PhononStructure::DisplacementPotential>();
	displPot_->build(
			this->get_primitive_unitcell_obj(),
			this->get_supercell_obj(),
			this->get_irrd_displmts_obj(),
			std::move(rsGridUC),
			std::move(rsGridSC),
			primitiveCellPotential,
			displPot);
}

void
ResourceHandler::initialize_wfct_obj()
{
	ElectronicStructure::Wavefunctions wfcts;
	std::string wfctdir = dataLoader_->get_optns().get_eld();
	if ( wfctdir.empty() )
		wfctdir = dataLoader_->get_optns().get_root_dir();
	wfcts.initialize(wfctdir, dataLoader_ );
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

} /* namespace IOMethods */
} /* namespace elephon */