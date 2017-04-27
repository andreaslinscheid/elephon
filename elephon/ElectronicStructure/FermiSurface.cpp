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

namespace elephon
{
namespace ElectronicStructure
{

FermiSurface::FermiSurface()
{

}

void FermiSurface::triangulate_Fermi_surface(
		std::vector<size_t> const& kgrid,
		size_t nbnd,
		std::vector<double> const& energies,
		size_t numTargetPoints)
{
	if ( kgrid.size() != 3 )
		throw std::logic_error("Can only triangulate 3D Fermi surfaces!");
	if ( kgrid[0]*kgrid[1]*kgrid[2]*nbnd != energies.size() )
		throw std::logic_error("Grid energy data size mismatch!");

	//Define the grid
	vtkSmartPointer<vtkStructuredPoints> grid =
			vtkSmartPointer<vtkStructuredPoints>::New();
	grid->SetOrigin(0,0,0);
	grid->SetDimensions(kgrid[0],kgrid[1],kgrid[2]);
	grid->SetSpacing(1.0/float(kgrid[0]),1.0/float(kgrid[1]),1.0/float(kgrid[2]));

	size_t dimGrid = kgrid[0]*kgrid[1]*kgrid[2];
	//Set data onto grid
	for ( size_t ib = 0 ; ib < nbnd; ib++)
	{
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

		// Decimation to reduce the number of triangles to roughly the numer set on input
		vtkSmartPointer<vtkDecimatePro> decimator =
				vtkDecimatePro::New();
		decimator->SetInputData(marched);
		decimator->SetTargetReduction(1.0-float(numTargetPoints)/float(marched->GetNumberOfPoints()));
		decimator->SetPreserveTopology(1);
		decimator->Update();

		decimator->GetOutput()->BuildCells();
		vtkSmartPointer<vtkIdList> pointIdsCell =
			  vtkSmartPointer<vtkIdList>::New();

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
		for ( size_t i = 0 ; i < decimator->GetOutput()->GetNumberOfPoints(); ++i )
			decimator->GetOutput()->GetPoint(i,&points[3*i]);

		kfPoints_.insert(std::end(kfPoints_), std::begin(points), std::end(points));
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

} /* namespace ElectronicStructure */
} /* namespace elephon */
