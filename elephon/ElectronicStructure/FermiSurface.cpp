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
#include <set>
#include <limits>

namespace elephon
{
namespace ElectronicStructure
{

void FermiSurface::triangulate_surface(
		LatticeStructure::RegularBareGrid grid,
		int nbnd,
		std::vector<double> const& energies,
		int numTargetPoints,
		double energyVal)
{
	grid_ = std::move(grid);
	auto d = grid_.get_grid_dim();
	if ( d[0]*d[1]*d[2]*nbnd != energies.size() )
		throw std::logic_error("Grid energy data size mismatch!");

	auto SurfaceGrid = d;
	for (auto &ki : SurfaceGrid )
		ki++;

	//Define the grid
	//We need to extend the grid to the first periodic point so that the surfaces will be correctly computed
	vtkSmartPointer<vtkImageData> dataGrid =
			vtkSmartPointer<vtkImageData>::New();
	dataGrid->SetOrigin(grid_.get_grid_shift()[0], grid_.get_grid_shift()[1], grid_.get_grid_shift()[2]);
	dataGrid->SetExtent(0,d[0],0,d[1],0,d[2]);
	dataGrid->SetSpacing(1.0/double(d[0]),1.0/double(d[1]),1.0/double(d[2]));

	//Loop the bands and compute the points, weights and gradient
	bandsMap_ = std::vector<int>( nbnd, -1 );
	vtkSmartPointer<vtkPolyData> * marched = new vtkSmartPointer<vtkPolyData> [nbnd] ;
	int totalNumberOfPointsFirstIteration = 0;
	for ( int ib = 0 ; ib < nbnd; ib++)
	{
		int surfaceGridDim = SurfaceGrid[0]*SurfaceGrid[1]*SurfaceGrid[2];
		//Set data onto grid. Use the periodicity for the last point in each dimension.
		vtkSmartPointer<vtkDoubleArray> doubleArray =
			  vtkSmartPointer<vtkDoubleArray>::New();
		doubleArray->SetNumberOfValues(surfaceGridDim);
		for (int k = 0 ; k < SurfaceGrid[2]; ++k)
			for (int j = 0 ; j < SurfaceGrid[1]; ++j)
				for (int i = 0 ; i < SurfaceGrid[0]; ++i)
				{
					//Note that vtk stores x,y,z in fast to slow running variables
					int consq = (k*SurfaceGrid[1]+j)*SurfaceGrid[0]+i;
					int ii = i%d[0];
					int jj = j%d[1];
					int kk = k%d[2];
					int consqEnergies = (kk*d[1]+jj)*d[0]+ii;
					doubleArray->SetValue(consq,energies[consqEnergies*nbnd+ib]);
				}
		dataGrid->GetPointData()->SetScalars( doubleArray );

		// Create a 3D model using marching cubes
		vtkSmartPointer<vtkMarchingCubes> mc =
				vtkSmartPointer<vtkMarchingCubes>::New();
		mc->SetInputData(dataGrid);
		mc->ComputeNormalsOff();
		mc->ComputeGradientsOff();//We need that for the adaptive weights
		mc->ComputeScalarsOff();
		// second value acts as threshold, and we are sampling the Fermi surface, i.e. E=0
		mc->SetValue(0, energyVal);

		// Create polydata from iso-surface
		mc->Update();
		marched[ib] = mc->GetOutput();
	}

	for ( int ib = 0 ; ib < nbnd; ib++)
		totalNumberOfPointsFirstIteration += marched[ib]->GetNumberOfPoints();

	double reductionRatio = 0.0;
	if ( totalNumberOfPointsFirstIteration > numTargetPoints )
	{
		reductionRatio = 1-double(numTargetPoints)/double(totalNumberOfPointsFirstIteration);
	}

	for ( int ib = 0 ; ib < nbnd; ib++)
	{
		if ( marched[ib]->GetNumberOfPoints() == 0 )
			continue;
		// Decimation to reduce the number of triangles to roughly the number set on input
		vtkSmartPointer<vtkDecimatePro> decimator =
				vtkSmartPointer<vtkDecimatePro>::New();
		decimator->SetInputData( marched[ib] );
		decimator->PreserveTopologyOff();
		decimator->BoundaryVertexDeletionOn();
		decimator->SetTargetReduction(reductionRatio);
		decimator->SetMaximumError( std::numeric_limits<double>::max() );
		decimator->SplittingOn();
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

			// in this additional step, we transport points to the actual isosurface according to a linear interpolation

			//transform to Cartesian coordinates
			grid_.get_lattice().reci_direct_to_cartesian_2pibya(p0,3);
			grid_.get_lattice().reci_direct_to_cartesian_2pibya(p1,3);
			grid_.get_lattice().reci_direct_to_cartesian_2pibya(p2,3);

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

int FermiSurface::get_npts_total() const
{
	return kfWeights_.size();
}

void FermiSurface::get_pt(int i, std::vector<double> & p) const
{
	if ( p.size() != 3 )
		p = std::vector<double>(3);
	std::copy(&kfPoints_[i*3],&kfPoints_[i*3]+3,p.data());
}

void
FermiSurface::get_pt_weight(int i, double & pw) const
{
	pw = kfWeights_[i];
}

int
FermiSurface::get_band_offset(int ib) const
{
	assert( (ib >= 0) && (ib < bandsMap_.size()) );
	return bandsMap_[ib];
}

std::vector<double> const&
FermiSurface::get_Fermi_vectors() const
{
	return kfPoints_;
}

std::vector<double> const&
FermiSurface::get_Fermi_weights() const
{
	return kfWeights_;
}

std::vector<double>
FermiSurface::get_Fermi_vectors_for_band(int ib) const
{
	int indexBandStart,indexBandEnd;
	this->band_index_range(ib,indexBandStart,indexBandEnd);
	std::vector<double> result( &kfPoints_[ indexBandStart*3 ], &kfPoints_[ indexBandEnd*3 ] );
	return result;
}

std::vector<double>
FermiSurface::get_Fermi_weights_for_band(int ib) const
{
	int indexBandStart,indexBandEnd;
	this->band_index_range(ib,indexBandStart,indexBandEnd);
	std::vector<double> result( &kfWeights_[ indexBandStart ], &kfWeights_[ indexBandEnd ] );
	return result;
}

void
FermiSurface::obtain_irreducible_Fermi_vectors_for_band(
		int ib,
		LatticeStructure::Symmetry const & symmetry,
		std::vector<double> & kpoints,
		std::vector<double> & weights) const
{
	int indexBandStart,indexBandEnd;
	this->band_index_range(ib,indexBandStart,indexBandEnd);
	assert(symmetry.is_reci());

	const double equalThresh = symmetry.get_symmetry_prec();
	typedef struct KPoint_{
		KPoint_(double x, double y, double z) : xi{x,y,z} {};
		std::vector<double> xi;
	} KPoint;

	auto compare_k_points = [&equalThresh] (KPoint const & k1, KPoint const & k2) {
		for ( int i = 0 ; i < 3; ++i)
			if ( std::abs(k1.xi[i] - k2.xi[i]) > equalThresh )
				return  k1.xi[i] < k2.xi[i];
		return false;
	};

	auto itKf = kfPoints_.begin();
	int nkf = indexBandEnd-indexBandStart;
	std::set<int> kfInIrred;
	for (int ikf = 0 ; ikf < nkf; ++ikf)
	{
		bool isIrreducibleZone = false;
		KPoint orig( *(itKf + ikf*3 ), *(itKf + ikf*3 + 1), *(itKf + ikf*3 + 2) );
		for ( int isym = 0 ; isym < symmetry.get_num_symmetries(); ++isym)
		{
			KPoint rot = orig;
			symmetry.rotate<double>(isym, rot.xi.begin(), rot.xi.end(), true);
			if ( compare_k_points(orig,rot) )
				break;
			// no rotated vector is smaller than the originial one
			isIrreducibleZone = true;
		}
		if ( isIrreducibleZone )
			kfInIrred.insert(ikf);
	}
	kpoints.resize( kfInIrred.size()*3 );
	weights.resize( kfInIrred.size() );

	auto it = kfInIrred.begin();
	for (int ikfirr = 0; ikfirr < kfInIrred.size(); ++ikfirr )
	{
		int ikf = *it++;
		std::copy(&kfPoints_[ikf*3], &kfPoints_[ikf*3]+3, &kpoints[ikfirr*3]);
		weights[ikfirr] = kfWeights_[ikf];
	}
}

void FermiSurface::band_index_range(
		int ib, int &start, int &end) const
{
	int numTotal = kfWeights_.size();
	start = end = 0;
	if ( bandsMap_[ib] >= 0 )
	{
		start = bandsMap_[ib];
		end = numTotal;
		if ( ib+1 < bandsMap_.size() )
			if (  bandsMap_[ib+1] > 0 )
				end = static_cast<int>(bandsMap_[ib+1]);
	}
}

} /* namespace ElectronicStructure */
} /* namespace elephon */
