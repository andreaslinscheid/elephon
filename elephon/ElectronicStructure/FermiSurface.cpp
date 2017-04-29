/*	This file FermiSurface.cpp is part of elephon.
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
 *  Created on: Apr 26, 2017
 *      Author: A. Linscheid
 */

#include <ElectronicStructure/FermiSurface.h>
#include <algorithm>
#include "vtkSmartPointer.h"
#include "vtkStructuredPoints.h"
#include "vtkMarchingCubes.h"
#include "vtkIdList.h"
#include "vtkTriangle.h"
#include "vtkPointData.h"
#include "vtkFloatArray.h"
#include "vtkDecimatePro.h"
#include "vtkMath.h"
#include <cmath>
#include <assert.h>
#include "Algorithms/TrilinearInterpolation.h"

namespace elephon
{
namespace ElectronicStructure
{

FermiSurface::FermiSurface()
{

}

void FermiSurface::triangulate_Fermi_surface(
		std::vector<size_t> kgrid,
		size_t nbnd,
		std::vector<double> const& energies,
		size_t numTargetPoints)
{
	kgrid_ = std::move(kgrid);
	if ( kgrid_.size() != 3 )
		throw std::logic_error("Can only triangulate 3D Fermi surfaces!");
	if ( kgrid_[0]*kgrid_[1]*kgrid_[2]*nbnd != energies.size() )
		throw std::logic_error("Grid energy data size mismatch!");

	//Define the grid
	vtkSmartPointer<vtkStructuredPoints> grid =
			vtkSmartPointer<vtkStructuredPoints>::New();
	grid->SetOrigin(0,0,0);
	grid->SetDimensions(kgrid_[0],kgrid_[1],kgrid_[2]);
	grid->SetSpacing(1.0/float(kgrid_[0]),1.0/float(kgrid_[1]),1.0/float(kgrid_[2]));

	//Define a kpoint with implicit ordering to construct a fast index table
	struct kpoint
	{
		kpoint(double x, double y, double z,size_t i) : x_(x),y_(y),z_(z),index_(i) {};

		bool operator> (kpoint const & kp) const
		{
			const double accuracy = 10e-8;
			auto cmp = [&] (double x, double y) { return std::fabs(x - y) > accuracy; };
			if ( cmp(x_,kp.x_) )
				return x_ < kp.x_;
			if ( cmp(y_,kp.y_) )
				return y_ < kp.y_;
			if ( cmp(z_,kp.z_) )
				return z_ < kp.z_;
			return false;
		}

		double x_,y_,z_;
		size_t index_;
	};

	//Loop the bands and compute the points, weights and gradient
	size_t dimGrid = kgrid_[0]*kgrid_[1]*kgrid_[2];
	bandsMap_ = std::vector<int>( nbnd, -1 );
	for ( size_t ib = 0 ; ib < nbnd; ib++)
	{
		//Set data onto grid
		vtkSmartPointer<vtkFloatArray> floatArray =
			  vtkSmartPointer<vtkFloatArray>::New();
		floatArray->SetNumberOfValues(energies.size());
		for (size_t i = 0 ; i < dimGrid; ++i)
			floatArray->InsertValue(i,energies[ib*dimGrid + i]);
		grid->GetPointData()->SetScalars( floatArray );

		// Create a 3D model using marching cubes
		vtkSmartPointer<vtkMarchingCubes> mc =
				vtkSmartPointer<vtkMarchingCubes>::New();
		mc->SetInputData(grid);
		mc->ComputeNormalsOn();
		mc->ComputeGradientsOn();
		// second value acts as threshold, and we are sampling the Fermi surface, i.e. E=0
		mc->SetValue(0, 0.0);

		// Create polydata from iso-surface
		vtkSmartPointer<vtkPolyData> marched =
				vtkSmartPointer<vtkPolyData>::New();
		mc->Update();
		marched->DeepCopy(mc->GetOutput());

		// Decimation to reduce the number of triangles to roughly the number set on input
		vtkSmartPointer<vtkDecimatePro> decimator =
				vtkDecimatePro::New();
		decimator->SetInputData(marched);
		decimator->SetTargetReduction(1.0-float(numTargetPoints)/float(marched->GetNumberOfPoints()));
		decimator->SetPreserveTopology(1);
		decimator->Update();

		decimator->GetOutput()->BuildCells();
		vtkSmartPointer<vtkIdList> pointIdsCell =
			  vtkSmartPointer<vtkIdList>::New();

		//compute the area of each rectangle and attribute 1/3 to each corner point.
		//Thus each point accumulates weight from rectangles it is part of.
		std::vector<double> kfWeightsBand(decimator->GetOutput()->GetNumberOfPoints(), 0.0);
		for(vtkIdType i = 0; i < decimator->GetOutput()->GetNumberOfCells(); i++)
		{
			vtkCell* cell = decimator->GetOutput()->GetCell(i);

			vtkTriangle* triangle = dynamic_cast<vtkTriangle*>(cell);
			double p0[3],p1[3],p2[3];
			triangle->GetPoints()->GetPoint(0, p0);
			triangle->GetPoints()->GetPoint(1, p1);
			triangle->GetPoints()->GetPoint(2, p2);

			double area = std::abs(vtkTriangle::TriangleArea(p0, p1, p2));

			decimator->GetOutput()->GetCellPoints(i,pointIdsCell);
			if ( pointIdsCell->GetNumberOfIds() != 3 )
				throw std::logic_error("Error, number of points not 3!");
			for (size_t j=0; j < 3; j++)
				kfWeightsBand[ pointIdsCell->GetId(j) ] += area/3.0;
		}
		kfWeights_.insert(std::end(kfWeights_), std::begin(kfWeightsBand), std::end(kfWeightsBand));

		std::vector<double> points(decimator->GetOutput()->GetNumberOfPoints()*3);
		for ( int i = 0 ; i < decimator->GetOutput()->GetNumberOfPoints(); ++i )
			decimator->GetOutput()->GetPoint(i,&points[3*i]);

		kfPoints_.insert(std::end(kfPoints_), std::begin(points), std::end(points));
		bandsMap_[ib] = kfWeights_.size();
	}

}

size_t FermiSurface::get_npts_total() const
{
	return kfWeights_.size();
}

void FermiSurface::get_pt(size_t i, std::vector<double> & p) const
{
	if ( p.size() != 3 )
		p = std::vector<double>(3);
	std::copy(&kfPoints_[i*3],&kfPoints_[i*3]+3,p.data());
}

void FermiSurface::get_pt_weight(size_t i, double & pw) const
{
	pw = kfWeights_[i];
}

void FermiSurface::compute_fermi_velocities(
		std::vector<double> const& energyGradientField)
{
	assert( energyGradientField.size() == bandsMap_.size()*kgrid_[0]*kgrid_[1]*kgrid_[2] );
	size_t nbnd = bandsMap_.size();

	//fetch the grid points that are required for a linear interpolation
	elephon::Algorithms::TrilinearInterpolation triLin( kgrid_ );
	std::vector<size_t> queryIndices;
	triLin.data_query( kfPoints_, queryIndices );

	//copy the relevant data
	std::vector<double> requestedData(queryIndices.size()*nbnd) ;
	auto itB = energyGradientField.begin();
	for (size_t i = 0 ; i < queryIndices.size(); ++i )
		std::copy(itB+queryIndices[i]*nbnd,itB+(queryIndices[i]+1)*nbnd,
				requestedData.begin()+i*nbnd);

	//interpolate
	triLin.interpolate( nbnd, requestedData, kfVeloc_);
}

} /* namespace ElectronicStructure */
} /* namespace elephon */
