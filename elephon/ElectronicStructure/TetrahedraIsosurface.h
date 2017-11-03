/*	This file TetrahedraIsosurface.h is part of elephon.
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
 *  Created on: Oct 9, 2017
 *      Author: A. Linscheid
 */

#ifndef ELEPHON_ELECTRONICSTRUCTURE_TETRAHEDRAISOSURFACE_H_
#define ELEPHON_ELECTRONICSTRUCTURE_TETRAHEDRAISOSURFACE_H_

#include "LatticeStructure/TetrahedraGrid.h"
#include "LatticeStructure/Tetrahedron.h"
#include "LatticeStructure/DataRegularGrid.h"
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <memory>
#include <vector>
#include <string>

namespace elephon
{
namespace ElectronicStructure
{

namespace detail {
struct NonGridPoint;
struct Triangle;
struct kinfo; };

class TetrahedraIsosurface
{
public:

	typedef detail::NonGridPoint NonGridPoint;

	typedef detail::Triangle Triangle;

	void initialize(
			std::shared_ptr<const LatticeStructure::TetrahedraGrid> tetra,
			std::shared_ptr<const LatticeStructure::DataRegularGrid<double>> bands,
			std::vector<double> const & isoEnergies);

	void write_polydata_file( std::string const & filename,
			int isoNum,
			std::shared_ptr<std::vector<double>> data = nullptr) const;

	void generate_vtk_polydata(
			int isoNum,
			vtkSmartPointer<vtkPolyData> & polydata) const;

	void get_irreducible_iso_vector_integration_weights(
			int isoNum,
			int ibnd,
			std::vector<double> & kIso,
			std::vector<double> & kIsoWeights) const;

	void get_irreducible_iso_vector_integration_weights_no_multiplicty(
			int isoNum,
			int ibnd,
			std::vector<double> & kIso,
			std::vector<double> & kIsoWeights) const;

	/**
	 *
	 *  Description of the Algorithm:
	 *
	 *
	 * @param isoNum
	 * @param ibnd
	 * @param kIso
	 * @param kIsoWeights
	 */
	void get_reducible_iso_vector_integration_weights(
			int isoNum,
			int ibnd,
			std::vector<double> & kIso,
			std::vector<double> & kIsoWeights) const;

	int get_num_iso_pts(int isoNum, int ibnd) const;

	int get_num_triangles(int isoNum, int ibnd) const;

	std::vector<Triangle> const & get_triangles(int isoNum, int ibnd) const;

	int get_nBnd() const;

	int get_num_iso_energies() const;

	std::shared_ptr<const LatticeStructure::TetrahedraGrid> get_tetra_grid() const;

private:

	int numIsoE_ = 0;

	int numBands_ = 0;

	std::vector<std::shared_ptr<std::vector<double>>> kIso_;

	/// dF / |grad E| surface integral weight in the irreducible cell
	/// This data array factors in the multiplicity of any given tetrahedron
	std::vector<std::vector<double>> kIsoWeights_;

	/// dF / |grad E| surface integral weight. This factors out the multiplicity
	std::vector<std::vector<double>> kIsoReducibleWeights_;

	std::vector<std::vector<int>> kIsoMultiplicity_;

	std::vector<std::vector<Triangle>> triangles_;

	std::shared_ptr<const LatticeStructure::TetrahedraGrid> tetra_;

	void compute_iso_vectors(
			LatticeStructure::Tetrahedron const & t,
			std::multimap<double,int> const & ecorner,
			double e,
			std::vector<double> & kIso) const;

	void set_multiplitcity(std::map<NonGridPoint, detail::kinfo> & irredIsoSurface);
};

namespace detail
{
struct NonGridPoint
{
	NonGridPoint(double kx, double ky, double kz, double gridprec)
		: k_{kx, ky, kz}, eps_(gridprec) { };

	const double k_[3];

	const double eps_;

	friend bool operator< (NonGridPoint const & ngp1, NonGridPoint const & ngp2);
};

struct Triangle
{
	Triangle(
			int ip1, int ip2, int ip3,
			std::shared_ptr<const std::vector<double>> points)
				: ip1_(ip1), ip2_(ip2), ip3_(ip3), pts_(points) { 	};

	const int ip1_, ip2_, ip3_;

	std::shared_ptr<const std::vector<double>> pts_;

	double area() const;

	std::vector<double> get_point_coords(int i) const;
};
};/* namespace detail */

} /* namespace ElectronicStructure */
} /* namespace elephon */

#endif /* ELEPHON_ELECTRONICSTRUCTURE_TETRAHEDRAISOSURFACE_H_ */
