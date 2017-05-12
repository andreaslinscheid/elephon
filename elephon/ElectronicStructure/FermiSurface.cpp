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
#include "vtkImageData.h"
#include "vtkMarchingCubes.h"
#include "vtkIdList.h"
#include "vtkTriangle.h"
#include "vtkPointData.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkDecimatePro.h"
#include "vtkMath.h"
#include <cmath>
#include <assert.h>

namespace elephon
{
namespace ElectronicStructure
{

FermiSurface::FermiSurface()
{

}

void FermiSurface::triangulate_surface(
		std::vector<size_t> kgrid,
		std::vector<double> const& reciprocalLatticeMatrix,
		size_t nbnd,
		std::vector<double> const& energies,
		size_t numTargetPoints,
		double energyVal)
{
	assert(reciprocalLatticeMatrix.size() == 9 );
	std::vector<double> const& B = reciprocalLatticeMatrix;
	auto to_cart = [&] (double pt[3])
	{
		double ptmp[3];
		std::copy(pt,pt+3,ptmp);
		for ( size_t i = 0; i < 3 ; ++i)
			pt[i] = (B[i*3+0]*ptmp[0]+B[i*3+1]*ptmp[1]+B[i*3+2]*ptmp[2]);
	};

	kgrid_ = std::move(kgrid);
	if ( kgrid_.size() != 3 )
		throw std::logic_error("Can only triangulate 3D Fermi surfaces!");
	if ( kgrid_[0]*kgrid_[1]*kgrid_[2]*nbnd != energies.size() )
		throw std::logic_error("Grid energy data size mismatch!");

	auto SurfaceGrid = kgrid_;
	for (auto &ki : SurfaceGrid )
		ki++;

	//Define the grid
	//We need to extend the grid to the first periodic point so that the surfaces will be correctly computed
	vtkSmartPointer<vtkImageData> grid =
			vtkSmartPointer<vtkImageData>::New();
	grid->SetOrigin(0,0,0);
	grid->SetExtent(0,kgrid_[0],0,kgrid_[1],0,kgrid_[2]);
	grid->SetSpacing(1.0/double(kgrid_[0]),1.0/double(kgrid_[1]),1.0/double(kgrid_[2]));

	//Loop the bands and compute the points, weights and gradient
	bandsMap_ = std::vector<int>( nbnd, -1 );
	vtkSmartPointer<vtkPolyData> * marched = new vtkSmartPointer<vtkPolyData> [nbnd] ;
	size_t totalNumberOfPointsFirstIteration = 0;
	for ( size_t ib = 0 ; ib < nbnd; ib++)
	{
		size_t surfaceGridDim = SurfaceGrid[0]*SurfaceGrid[1]*SurfaceGrid[2];
		//Set data onto grid. Use the periodicity for the last point in each dimension.
		vtkSmartPointer<vtkDoubleArray> doubleArray =
			  vtkSmartPointer<vtkDoubleArray>::New();
		doubleArray->SetNumberOfValues(surfaceGridDim);
		for (size_t i = 0 ; i < SurfaceGrid[0]; ++i)
			for (size_t j = 0 ; j < SurfaceGrid[1]; ++j)
				for (size_t k = 0 ; k < SurfaceGrid[2]; ++k)
				{
					//Note that vtk stores x,y,z in fast to slow running variables
					size_t consq = (k*SurfaceGrid[1]+j)*SurfaceGrid[0]+i;
					size_t ii = i%kgrid_[0];
					size_t jj = j%kgrid_[1];
					size_t kk = k%kgrid_[2];
					size_t consqEnergies = (ii*kgrid_[1]+jj)*kgrid_[2]+kk;
					doubleArray->SetValue(consq,energies[consqEnergies*nbnd+ib]);
				}
		grid->GetPointData()->SetScalars( doubleArray );

		// Create a 3D model using marching cubes
		vtkSmartPointer<vtkMarchingCubes> mc =
				vtkSmartPointer<vtkMarchingCubes>::New();
		mc->SetInputData(grid);
		mc->ComputeNormalsOff();
		mc->ComputeGradientsOff();//We need that for the adaptive weights
		mc->ComputeScalarsOff();
		// second value acts as threshold, and we are sampling the Fermi surface, i.e. E=0
		mc->SetValue(0, energyVal);

		// Create polydata from iso-surface
		mc->Update();
		marched[ib] = mc->GetOutput();
	}

	for ( size_t ib = 0 ; ib < nbnd; ib++)
		totalNumberOfPointsFirstIteration += marched[ib]->GetNumberOfPoints();

	double reductionRatio = 1.0-double(numTargetPoints)/double(totalNumberOfPointsFirstIteration);
	reductionRatio -= std::floor(reductionRatio);

	for ( size_t ib = 0 ; ib < nbnd; ib++)
	{
		if ( marched[ib]->GetNumberOfPoints() == 0 )
			continue;
		// Decimation to reduce the number of triangles to roughly the number set on input
		vtkSmartPointer<vtkDecimatePro> decimator =
				vtkSmartPointer<vtkDecimatePro>::New();
		decimator->SetInputData( marched[ib] );
		decimator->SetTargetReduction(reductionRatio);
		decimator->SetPreserveTopology(1);
		decimator->Update();

		decimator->GetOutput()->BuildCells();
		vtkSmartPointer<vtkIdList> pointIdsCell =
			  vtkSmartPointer<vtkIdList>::New();

		int npts = decimator->GetOutput()->GetNumberOfPoints();
		int ncells = decimator->GetOutput()->GetNumberOfCells();

		std::vector<double> kfWeightsBand(npts, 0.0);
		for(vtkIdType i = 0; i < ncells; i++)
		{
			double p0[3],p1[3],p2[3],center[3],p01[3],p02[3],p12[3];
			decimator->GetOutput()->GetCellPoints(i,pointIdsCell);
			if ( pointIdsCell->GetNumberOfIds() != 3 )
				throw std::logic_error("Error, number of points not 3!");

			decimator->GetOutput()->GetPoint(pointIdsCell->GetId(0),p0);
			decimator->GetOutput()->GetPoint(pointIdsCell->GetId(1),p1);
			decimator->GetOutput()->GetPoint(pointIdsCell->GetId(2),p2);

			//transform to Cartesian coordinates
			to_cart(p0);
			to_cart(p1);
			to_cart(p2);

			for (int xi = 0 ; xi < 3; xi++)
			{
				center[xi] = 1.0/3.0*(p0[xi]+p1[xi]+p2[xi]);
				p01[xi] = 1.0/2.0*(p0[xi]+p1[xi]);
				p02[xi] = 1.0/2.0*(p0[xi]+p2[xi]);
				p12[xi] = 1.0/2.0*(p1[xi]+p2[xi]);
			}

			//Get triangles area
			auto dotProd = [] (const double p1[3], const double p2[3])
			{
					return (p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1])
							+(p1[2]-p2[2])*(p1[2]-p2[2]);
			};
			auto areaTri = [&] (const double p1[3], const double p2[3], const double p3[3])
			{
				double a = dotProd(p1,p2);
				double b = dotProd(p2,p3);
				double c = dotProd(p3,p1);
				return (0.25* std::sqrt(std::fabs(4.0*a*c - (a-b+c)*(a-b+c))));
			};

			//Attribute weights to points
			kfWeightsBand[ pointIdsCell->GetId(0) ] += areaTri(p0,p01,center)
													+ areaTri(p0,p02,center);
			kfWeightsBand[ pointIdsCell->GetId(1) ] += areaTri(p1,p01,center)
													+ areaTri(p1,p12,center);
			kfWeightsBand[ pointIdsCell->GetId(2) ] += areaTri(p2,p02,center)
													+ areaTri(p2,p12,center);
		}
		bandsMap_[ib] = kfWeights_.size();
		kfWeights_.insert(std::end(kfWeights_),std::begin(kfWeightsBand),std::end(kfWeightsBand));

		std::vector<double> points(decimator->GetOutput()->GetNumberOfPoints()*3);
		for ( int ip = 0 ; ip < decimator->GetOutput()->GetNumberOfPoints(); ++ip )
			decimator->GetOutput()->GetPoint(ip,&points[3*ip]);

		kfPoints_.insert(std::end(kfPoints_),std::begin(points),std::end(points));
	}
	delete [] marched;
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

std::vector<double> const&
FermiSurface::get_Fermi_vectors() const
{
	return kfPoints_;
}

std::vector<double>
FermiSurface::get_Fermi_vectors_for_band(size_t ib) const
{
	size_t indexBandStart,indexBandEnd;
	this->band_index_range(ib,indexBandStart,indexBandEnd);
	std::vector<double> result( &kfPoints_[ indexBandStart*3 ], &kfPoints_[ indexBandEnd*3 ] );
	return result;
}

std::vector<double>
FermiSurface::get_Fermi_weights_for_band(size_t ib) const
{
	size_t indexBandStart,indexBandEnd;
	this->band_index_range(ib,indexBandStart,indexBandEnd);
	std::vector<double> result( &kfWeights_[ indexBandStart ], &kfWeights_[ indexBandEnd ] );
	return result;
}

void FermiSurface::band_index_range(
		size_t ib, size_t &start, size_t &end) const
{
	size_t numTotal = kfWeights_.size();
	start = end = 0;
	if ( bandsMap_[ib] >= 0 )
	{
		start = bandsMap_[ib];
		end = numTotal;
		if ( ib+1 < bandsMap_.size() )
			if (  bandsMap_[ib+1] > 0 )
				end = static_cast<size_t>(bandsMap_[ib+1]);
	}
}

} /* namespace ElectronicStructure */
} /* namespace elephon */
