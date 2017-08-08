/*	This file LatticeModule.h is part of elephon.
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
 *  Created on: May 15, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_LATTICESTRUCTURE_LATTICEMODULE_H_
#define ELEPHON_LATTICESTRUCTURE_LATTICEMODULE_H_

#include <vector>

namespace elephon
{
namespace LatticeStructure
{

class LatticeModule
{
public:

	LatticeModule();

	LatticeModule( std::vector<double> latticeMatrix );

	/**
	 * latticeMatrix : Cartesian, units are in Angstrom
	 */
	void initialize( std::vector<double> latticeMatrix );

	std::vector<double> const & get_latticeMatrix() const;

	std::vector<double> const & get_reciprocal_latticeMatrix() const;

	double get_volume() const;

	double get_alat() const;

	std::vector<double> get_lattice_vector(int n) const;

	std::vector<double> get_reci_lattice_vector(int n) const;

	void direct_to_cartesian_matrix(double * mat, int nelem) const;

	void reci_cartesian_to_direct_matrix(double * mat, int nelem) const;

	void direct_to_cartesian(double p[3]) const;

	void direct_to_cartesian(double * p, int nelem) const;

	void direct_to_cartesian(std::vector<double> & v) const;

	void direct_to_cartesian_angstroem(std::vector<double> & v) const;

	void direct_to_cartesian_angstroem(double * mat, int nelem) const;

	void cartesian_to_direct(std::vector<double> & v) const;

	void reci_direct_to_cartesian(std::vector<double> & v) const;

	void reci_direct_to_cartesian(double * p , int nelem) const;

	void reci_direct_to_cartesian_2pibya(double * p , int nelem) const;

	void reci_direct_to_cartesian_2pibya(std::vector<double> & v) const;
private:

	double vol_ = 0;

	double alat_ = 0;

	//units are in alat_
	std::vector<double> latticeMatrix_;

	//units are in 1/alat_
	std::vector<double> reciLatMatrix_;

	void cross_prod( std::vector<double> const& v1,
			std::vector<double> const& v2,
			std::vector<double> & v1xv2) const;
};

} /* namespace LatticeStructure */
} /* namespace elephon */

#endif /* ELEPHON_LATTICESTRUCTURE_LATTICEMODULE_H_ */
