/*	This file AtomDisplacementCollection.cpp is part of elephon.
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
 *  Created on: Feb 10, 2018
 *      Author: A. Linscheid
 */

#include "LatticeStructure/AtomDisplacementCollection.h"
#include "LatticeStructure/AtomSymmetryConnection.h"
#include "LatticeStructure/AtomDisplacement.h"
#include "LatticeStructure/UnitCell.h"
#include "LatticeStructure/AtomSymmetryConnection.h"

namespace elephon
{
namespace LatticeStructure
{

void
AtomDisplacementCollection::initialize(
		std::shared_ptr<const UnitCell> primitiveCell,
		bool symmetricDisplacements,
		double displacementMagnitude)
{
	symmetricDisplacements_ = symmetricDisplacements;
	displMagn_ = displacementMagnitude;
	atomSym_ = primitiveCell->get_atom_symmetry();
	this->generate_displacements(primitiveCell, symmetricDisplacements_, displMagn_);
}

std::string
AtomDisplacementCollection::get_relative_folder_structure_displ_run(int totalIrredIndex) const
{
	assert((totalIrredIndex>=0)&&(totalIrredIndex<this->get_tota_num_irred_displacements()));
	auto p = get_total_irred_index_to_atom_and_rel_irred_index(totalIrredIndex);
	return std::string("displ_")+std::to_string(p.first)+"_"+std::to_string(p.second);
}

std::pair<int,int>
AtomDisplacementCollection::get_total_irred_index_to_atom_and_rel_irred_index(int totalIrredIndex) const
{
	assert((totalIrredIndex>=0) && (totalIrredIndex<totalIrredToAtomAndRelIrred_.size()));
	return totalIrredToAtomAndRelIrred_[totalIrredIndex];
}

int
AtomDisplacementCollection::get_tota_num_irred_displacements() const
{
	int sum = 0;
	for (int ia = 0; ia < irredDispl_.size(); ++ia)
		sum += irredDispl_[ia].size();
	return sum;
}

bool
AtomDisplacementCollection::check_atom_is_irreducible(int atomIndex) const
{
	return atomSym_->check_atom_is_irreducible(atomIndex);
}

int
AtomDisplacementCollection::get_num_irred_displacements_for_atom(int atomIndex) const
{
	assert((atomIndex>=0)&&(atomIndex<mapReducibleToIrreducibleAtoms_.size()));
	const int irredIndex = mapReducibleToIrreducibleAtoms_[atomIndex];
	assert((irredIndex>=0)&&(irredIndex<irredDispl_.size()));
	return irredDispl_[irredIndex].size();
}

int
AtomDisplacementCollection::get_num_atoms_primitive_cell() const
{
	return mapReducibleToIrreducibleAtoms_.size();
}

std::vector<AtomDisplacement>
AtomDisplacementCollection::get_red_displacements_for_atom(int atomIndex) const
{
	assert((atomIndex>=0)&&(atomIndex<mapReducibleToIrreducibleAtoms_.size()));
	const int irredIndex = mapReducibleToIrreducibleAtoms_[atomIndex];
	assert((irredIndex>=0)&&(irredIndex<irredDispl_.size()));
	return redDispl_[irredIndex];
}

std::vector<std::pair<int, std::vector<AtomDisplacement>>>
AtomDisplacementCollection::get_irreducible_displacements() const
{
	std::vector<std::pair<int, std::vector<AtomDisplacement>>> result(irredDispl_.size());
	for (int ia = 0; ia < irredDispl_.size(); ++ia)
	{
		result[ia].first = atomSym_->get_list_irreducible_atoms()[ia];
		result[ia].second = irredDispl_[ia];
	}
	return result;
}

int
AtomDisplacementCollection::get_num_red_displacements_for_atom(int atomIndex) const
{
	assert((atomIndex>=0)&&(atomIndex<redDispl_.size()));
	return redDispl_[atomIndex].size();
}

std::vector<std::pair<int,int>>
AtomDisplacementCollection::get_symmetry_relation_red_displacements_for_atom(int atomIndex) const
{
	assert((atomIndex>=0) && (atomIndex<mapReducibleToIrreducibleAtoms_.size()));
	const int irredAtomIndex = mapReducibleToIrreducibleAtoms_[atomIndex];
	const int nRedDispl = redDispl_[irredAtomIndex].size();
	std::vector<std::pair<int,int>> result(nRedDispl);
	for (int ird = 0; ird < nRedDispl; ++ird)
	{
		result[ird].first = redToIrredDispl_[irredAtomIndex][ird];
		const int symRedToIrred = symRedToIrredDispl_[irredAtomIndex][ird];
		result[ird].second = symIrredToRedDispl_[irredAtomIndex][result[ird].first][symRedToIrred];
	}
	return result;

}

void
AtomDisplacementCollection::generate_pseudo_inverse_reducible_displacements_for_atom(
		int atomIndex,
		Auxillary::Multi_array<double,2> & pInv) const
{
	//build the transpose of the U matrix from the irreducible displacements
	auto reducibleDisplacements = this->get_red_displacements_for_atom(atomIndex);
	const int nRedDispl = reducibleDisplacements.size();
	Auxillary::alignedvector::DV U( 3 * nRedDispl ), pInvU;
	for ( int idr = 0 ; idr < nRedDispl;  ++idr)
	{
		auto d = reducibleDisplacements[idr].get_direction();
		for ( auto &xi : d )
			xi *= reducibleDisplacements[idr].get_magnitude();
		std::copy( d.begin(), d.end(),  &U[idr*3] );
	}

	//compute the pseudo inverse
	Algorithms::LinearAlgebraInterface linAlg;
	linAlg.pseudo_inverse( std::move(U), nRedDispl, 3, pInvU );
	pInv.resize(boost::extents[3][nRedDispl]);
	std::copy(pInvU.begin(), pInvU.end(), pInv.data());
}

bool
AtomDisplacementCollection::get_treat_displacements_pm_symmetric() const
{
	return symmetricDisplacements_;
}

void
AtomDisplacementCollection::generate_displacements(
		std::shared_ptr<const UnitCell> primitiveCell,
		bool symmetricDisplacements,
		double displMagn)
{
	irredDispl_.clear();
	redDispl_.clear();
	symIrredToRedDispl_.clear();
	irredToRedDispl_.clear();
	symRedToIrredDispl_.clear();
	redToIrredDispl_.clear();

	const int nIrredAtoms = atomSym_->get_list_irreducible_atoms().size();
	irredDispl_.reserve(nIrredAtoms);
	redDispl_.reserve(nIrredAtoms);
	symIrredToRedDispl_.reserve(nIrredAtoms);
	irredToRedDispl_.reserve(nIrredAtoms);
	symRedToIrredDispl_.reserve(nIrredAtoms);
	redToIrredDispl_.reserve(nIrredAtoms);

	int irreducibleIndex = 0;
	totalIrredToAtomAndRelIrred_.clear();
	mapReducibleToIrreducibleAtoms_.resize(primitiveCell->get_atoms_list().size());
	for ( int atomIndex : atomSym_->get_list_irreducible_atoms())
	{
		// mark all atom indices in the star as pointing to this irreducible index
		for (auto star : atomSym_->get_star_atom_indices(atomIndex))
			mapReducibleToIrreducibleAtoms_[star.first] = irreducibleIndex;

		std::vector<AtomDisplacement> irreducibleThisAtom, reducibleThisAtom;
		std::vector<int> redToIrred, symRedToIrred;
		std::vector< std::vector<int> > irredToRed, symIrredToRed;

		this->get_site_displacements(
				primitiveCell->get_lattice(),
				primitiveCell->get_atoms_list()[atomIndex],
				symmetricDisplacements, primitiveCell->get_site_symmetry(atomIndex), displMagn,
				irreducibleThisAtom, reducibleThisAtom,
				 redToIrred, symRedToIrred,
				 irredToRed, symIrredToRed);

		for ( int iRelDispl = 0 ; iRelDispl < irreducibleThisAtom.size(); ++iRelDispl)
			totalIrredToAtomAndRelIrred_.push_back(std::make_pair(atomIndex, iRelDispl));

		irredDispl_.push_back(std::move(irreducibleThisAtom));
		redDispl_.push_back(std::move(reducibleThisAtom));
		symIrredToRedDispl_.push_back(std::move(symIrredToRed));
		irredToRedDispl_.push_back(std::move(irredToRed));
		symRedToIrredDispl_.push_back(std::move(symRedToIrred));
		redToIrredDispl_.push_back(std::move(redToIrred));
		irreducibleIndex++;
	}
}

void
AtomDisplacementCollection::get_site_displacements(
		LatticeStructure::LatticeModule const & lattice,
		Atom const & atomicSite,
		bool symmetricDisplacements,
		LatticeStructure::Symmetry const & siteSymmetry,
		double displMagn,
		std::vector<AtomDisplacement> & irreducible,
		std::vector<AtomDisplacement> & reducible,
		std::vector<int> & redToIrredDispl,
		std::vector<int> & symRedToIrredDispl,
		std::vector< std::vector<int> > & irredToRedDispl,
		std::vector< std::vector<int> > & symIrredToRedDispl) const
{
	typedef std::vector<double> V;
	auto pos = atomicSite.get_position();
	auto name = atomicSite.get_kind();

	double delta = siteSymmetry.get_symmetry_prec();
	Algorithms::LinearAlgebraInterface linAlg;

	//The best strategy to save numerical effort is to reduce the number of irreducible displacements.
	//Thus, one should take displacements that are mapped to linear independent displacements by a site symmetry.
	//On the other hand, in order to make the calculation of the given displacement more efficient, the
	//overall symmetry group should not be affected.
	std::vector<AtomDisplacement> nonRotated;
	std::vector<AtomDisplacement> rotated;

	auto symmetry_expand = [&] () {
		rotated.clear();
		rotated.reserve( nonRotated.size() * siteSymmetry.get_num_symmetries() );
		irredToRedDispl.resize(nonRotated.size(), std::vector<int>( siteSymmetry.get_num_symmetries() ));
		symIrredToRedDispl.resize(nonRotated.size(), std::vector<int>( siteSymmetry.get_num_symmetries() ));
		redToIrredDispl.resize(nonRotated.size()* siteSymmetry.get_num_symmetries());
		symRedToIrredDispl.resize(nonRotated.size()* siteSymmetry.get_num_symmetries());
		for ( int issym = 0 ; issym < siteSymmetry.get_num_symmetries(); ++issym )
		{
			for (int inr = 0 ; inr < nonRotated.size(); ++inr)
			{
				irredToRedDispl[inr][issym] = rotated.size();
				symIrredToRedDispl[inr][issym] = issym;
				redToIrredDispl[rotated.size()] = inr;
				symRedToIrredDispl[rotated.size()] = siteSymmetry.get_index_inverse(issym);
				auto ad = nonRotated[inr];
				ad.transform_direction( siteSymmetry.get_sym_op(issym) );
				rotated.push_back( std::move(ad) );
			}
		}
	};

	//Insert displacement along lattice x
	V vx = {1,0,0};
	lattice.direct_to_cartesian_angstroem(vx);
	nonRotated.push_back( AtomDisplacement(name, displMagn, pos, vx, delta, (not symmetricDisplacements) ) );

	//expand to all reducible displacements at this site
	symmetry_expand();

	//Add the minus displacement if it is not already in the set
	V mvx = {-1,0,0};
	lattice.direct_to_cartesian_angstroem(mvx);
	AtomDisplacement minusDx(name, displMagn, pos, mvx, delta, (not symmetricDisplacements) );
	auto it = std::find( rotated.begin(), rotated.end(), minusDx);
	if ( (it == rotated.end()) and (symmetricDisplacements) )
	{
		nonRotated.push_back( std::move(minusDx) );
		symmetry_expand();
	}

	V redDisplTmp;
	redDisplTmp.resize(3*rotated.size());
	for ( int i = 0 ; i < rotated.size(); ++i)
		for ( int j = 0 ; j < 3; ++j)
		redDisplTmp[i*3+j] = rotated[i].get_direction()[j];

	//Check if the resulting set of vectors spans full 3D space.
	//Very small overlaps (<0.01) are considered null
	int nullDim;
	V nullSpace;
	linAlg.null_space(redDisplTmp,rotated.size(),3,nullDim,nullSpace, 0.01);

	auto check_overlap_is_small = [] ( V const & nullSpace, int nullDim, V const & trialVector) {
		const double small = 0.01;
		assert(nullSpace.size() == 3*nullDim );
		assert(trialVector.size() == 3 );
		bool overlapIsZero = true;
		for ( int i = 0 ; i < nullDim ; ++i )
		{
			double overlap = 0;
			for ( int j = 0 ; j < 3 ; ++j )
				overlap += nullSpace[i*3+j]*trialVector[j];
			overlapIsZero = overlapIsZero and (overlap < small);
		}
		return overlapIsZero;
	};

	if ( nullDim != 0 )
	{
		//x + rotations do not span 3D space, insert displacement
		// along y (and -y) if its overlap with the nullSpace is not zero
		V vy = {0,1,0};
		lattice.direct_to_cartesian_angstroem(vy);
		if ( not check_overlap_is_small(nullSpace, nullDim, vy) )
			nonRotated.push_back( AtomDisplacement(name, displMagn, pos, vy, delta, (not symmetricDisplacements) ) );

		//nonRotated contains now x and y - expand by symmetry and check null space again
		symmetry_expand();

		V mvy = {0,-1,0};
		lattice.direct_to_cartesian_angstroem(mvy);
		AtomDisplacement minusDy(name, displMagn, pos, mvy, delta, (not symmetricDisplacements) );
		auto it = std::find( rotated.begin(), rotated.end(), minusDy);
		if ( (it == rotated.end()) and (symmetricDisplacements) )
		{
			nonRotated.push_back( std::move(minusDy) );
			symmetry_expand();
		}

		redDisplTmp.resize(3*rotated.size());
		for ( int i = 0 ; i < rotated.size(); ++i)
			for ( int j = 0 ; j < 3; ++j)
			redDisplTmp[i*3+j] = rotated[i].get_direction()[j];
		linAlg.null_space(redDisplTmp,rotated.size(),3,nullDim,nullSpace, 0.01);
	}

	if ( nullDim != 0 )
	{
		//x, y + rotations do not span 3D space, insert displacement
		// along z (and -z) if its overlap with the nullSpace is not zero
		V vz = {0,0,1};
		lattice.direct_to_cartesian_angstroem(vz);
		nonRotated.push_back( AtomDisplacement(name, displMagn, pos, vz, delta, (not symmetricDisplacements) ) );

		symmetry_expand();

		V mvz = {0,0,-1};
		lattice.direct_to_cartesian_angstroem(mvz);
		AtomDisplacement minusDz(name, displMagn, pos, mvz, delta, (not symmetricDisplacements) );
		auto it = std::find( rotated.begin(), rotated.end(), minusDz);
		if ( (it == rotated.end()) and (symmetricDisplacements) )
		{
			nonRotated.push_back( std::move(minusDz) );
			symmetry_expand();
		}

		//In principle, if the algorithm work, the following is unnecessary but better check
		//nonRotated contains now x, y and z - expand by symmetry and check null space again
		redDisplTmp.resize(3*rotated.size());
		for ( int i = 0 ; i < rotated.size(); ++i)
			for ( int j = 0 ; j < 3; ++j)
			redDisplTmp[i*3+j] = rotated[i].get_direction()[j];
		linAlg.null_space(redDisplTmp,rotated.size(),3,nullDim,nullSpace, 0.01);
	}

	if ( nullDim != 0 )
		throw std::logic_error( "Displacement generation algorithm failed to span all 3D space!" );

	irreducible  = std::move(nonRotated);
	reducible = std::move(rotated);
}

} /* namespace LatticeStructure */
} /* namespace elephon */
