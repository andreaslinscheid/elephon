/*	This file test_TetrahedraIsosurface.cpp is part of elephon.
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
 *  Created on: Oct 10, 2017
 *      Author: A. Linscheid
 */

#define BOOST_TEST_MODULE TetraGrids
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "ElectronicStructure/TetrahedraIsosurface.h"
#include "fixtures/DataLoader.h"
#include "fixtures/MockStartup.h"
#include <memory>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>

void plot_iso_surface(
		std::shared_ptr<const elephon::ElectronicStructure::TetrahedraIsosurface> iso,
		int ie,
		bool plotReducibleZone = false)
{
	// Create a polydata object and add everything to it
	vtkSmartPointer<vtkPolyData> polydata;
	iso->generate_vtk_polydata(ie, polydata);

	// Visualize
	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	mapper->SetInputConnection(polydata->GetProducerPort());
#else
	mapper->SetInputData(polydata);
#endif

	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);

	vtkSmartPointer<vtkRenderer> renderer =
	vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow =
	vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
	vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	renderer->AddActor(actor);

	renderWindow->Render();
	renderWindowInteractor->Start();
}

BOOST_AUTO_TEST_CASE( Sphere_simple )
{
	//construct data which will yield a sphere of radius r
	double const r = 0.41;
	auto grid = std::make_shared<elephon::LatticeStructure::RegularSymmetricGrid>();
	elephon::LatticeStructure::Symmetry sym;
	sym.set_reciprocal_space_sym();
	grid->initialize( {30,30,30}, 1e-6, {0.0, 0.0, 0.0}, sym, elephon::LatticeStructure::LatticeModule() );
	auto d = grid->get_grid_dim();
	elephon::Auxillary::alignedvector::DV data(d[0]*d[1]*d[2]);
	for (int k = 0 ; k < d[2]; ++k)
	  for (int j = 0 ; j < d[1]; ++j)
		for (int i = 0 ; i < d[0]; ++i)
		  {
			  double x = (i<=d[0]/2?double(i):double(i)-d[0])/double(d[0]);
			  double y = (j<=d[1]/2?double(j):double(j)-d[1])/double(d[1]);
			  double z = (k<=d[2]/2?double(k):double(k)-d[2])/double(d[2]);
			  data[(k*d[1]+j)*d[0]+i] = std::sqrt(x*x+y*y+z*z)/r;
		  }
	auto bands = std::make_shared<elephon::ElectronicStructure::ElectronicBands>();
	bands->initialize(1, 0.0, std::move(data), *grid);

	auto tetra = std::make_shared<elephon::LatticeStructure::TetrahedraGrid>();
	tetra->initialize(grid);

	auto tetraIso = std::make_shared<elephon::ElectronicStructure::TetrahedraIsosurface>();
	tetraIso->initialize(tetra, bands, {1.0});

	double surfaceArea = 0;
	for ( auto const & t : tetraIso->get_triangles(0, 0) )
		surfaceArea += t.area();

	std::cout << "Surface area of a sphere (triangulated): "
			<< surfaceArea << " expected: " << 4*M_PI*r*r << std::endl;
	BOOST_CHECK( std::fabs(surfaceArea-4*M_PI*r*r) < 1e-2);
	// uncomment for manual inspection
//	plot_iso_surface(tetraIso, 0);
}

BOOST_AUTO_TEST_CASE( wall_simple )
{
	auto grid = std::make_shared<elephon::LatticeStructure::RegularSymmetricGrid>();
	elephon::LatticeStructure::Symmetry sym;
	sym.set_reciprocal_space_sym();
	grid->initialize( {2, 1, 1}, 1e-6, {0.0, 0.0, 0.0}, sym, elephon::LatticeStructure::LatticeModule() );
	auto d = grid->get_grid_dim();
	elephon::Auxillary::alignedvector::DV data(d[0]*d[1]*d[2]);
	for (int k = 0 ; k < d[2]; ++k)
		for (int j = 0 ; j < d[1]; ++j)
			for (int i = 0 ; i < d[0]; ++i)
				data[(k*d[1]+j)*d[0]+i] = 2.0*(i==1?1:0);
	auto bands = std::make_shared<elephon::ElectronicStructure::ElectronicBands>();
	bands->initialize(1, 0.0, std::move(data), *grid);

	auto tetra = std::make_shared<elephon::LatticeStructure::TetrahedraGrid>();
	tetra->initialize(grid);

	auto tetraIso = std::make_shared<elephon::ElectronicStructure::TetrahedraIsosurface>();
	tetraIso->initialize(tetra, bands, {1.0});

	double surfaceArea = 0;
	for ( auto const & t : tetraIso->get_triangles(0, 0) )
		surfaceArea += t.area();

	std::cout << "Surface area of a double wall (triangulated): "
			<< surfaceArea << " expected: " << 2 << std::endl;
	BOOST_CHECK( std::fabs(surfaceArea- 2) < 1e-2);

	for ( auto const & t : tetraIso->get_triangles(0, 0) )
	{
		if ( t.get_point_coords(0)[0] < 0.5 )
		{
			BOOST_CHECK_CLOSE(t.get_point_coords(0)[0], 0.25, 0.1);
			BOOST_CHECK_CLOSE(t.get_point_coords(1)[0], 0.25, 0.1);
			BOOST_CHECK_CLOSE(t.get_point_coords(2)[0], 0.25, 0.1);
		}
		else
		{
			BOOST_CHECK_CLOSE(t.get_point_coords(0)[0], 0.75, 0.1);
			BOOST_CHECK_CLOSE(t.get_point_coords(1)[0], 0.75, 0.1);
			BOOST_CHECK_CLOSE(t.get_point_coords(2)[0], 0.75, 0.1);
		}
	}

	// uncomment for manual inspection
//	plot_iso_surface(tetraIso, 0);
}

BOOST_AUTO_TEST_CASE( box_simple )
{
	auto grid = std::make_shared<elephon::LatticeStructure::RegularSymmetricGrid>();
	elephon::LatticeStructure::Symmetry sym;
	sym.set_reciprocal_space_sym();
	grid->initialize( {4,4,4}, 1e-6, {0.0, 0.0, 0.0}, sym, elephon::LatticeStructure::LatticeModule() );
	auto d = grid->get_grid_dim();
	elephon::Auxillary::alignedvector::DV data(d[0]*d[1]*d[2]);
	for (int k = 0 ; k < d[2]; ++k)
		for (int j = 0 ; j < d[1]; ++j)
			for (int i = 0 ; i < d[0]; ++i)
			{
				  double x = (i<d[0]/2?double(i):double(i)-d[0])/double(d[0]);
				  double y = (j<d[1]/2?double(j):double(j)-d[1])/double(d[1]);
				  double z = (k<d[2]/2?double(k):double(k)-d[2])/double(d[2]);
				data[(k*d[1]+j)*d[0]+i] = 2.0*(fabs(x)>=0.25?1:0)*(fabs(y)>=0.25?1:0)*(fabs(z)>=0.25?1:0);
			}
	auto bands = std::make_shared<elephon::ElectronicStructure::ElectronicBands>();
	bands->initialize(1, 0.0, std::move(data), *grid);

	auto tetra = std::make_shared<elephon::LatticeStructure::TetrahedraGrid>();
	tetra->initialize(grid);

	auto tetraIso = std::make_shared<elephon::ElectronicStructure::TetrahedraIsosurface>();
	tetraIso->initialize(tetra, bands, {2.0});

	double surfaceArea = 0;
	for ( auto const & t : tetraIso->get_triangles(0, 0) )
		surfaceArea += t.area();

	std::cout << "Surface area of a box (triangulated): "
			<< surfaceArea << " expected: " << 6*0.25 << std::endl;
	BOOST_CHECK( std::fabs(surfaceArea- 6*0.25) < 1e-2);

	// uncomment for manual inspection
//	plot_iso_surface(tetraIso, 0);
}

BOOST_AUTO_TEST_CASE( TetrahedraIsosurface_Al )
{
	elephon::test::fixtures::MockStartup ms;
	auto rootDir = ms.get_data_for_testing_dir() / "Al" / "vasp" / "fcc_primitive";
	auto surfDataFile = rootDir / "surface.vtp";
	boost::filesystem::remove(surfDataFile);
	elephon::test::fixtures::DataLoader dl;
	auto resHndler = dl.create_resource_handler(std::string()+
			"root_dir = "+rootDir.string()
			);
	auto bands = resHndler->get_electronic_bands_obj();

	std::vector<int> bndIndices(bands->get_nBnd());
	for (int ibnd = 0 ; ibnd < bndIndices.size(); ++ibnd)
		bndIndices[ibnd] = ibnd;

	elephon::Auxillary::alignedvector::DV reducibleBndData;
	bands->generate_reducible_data(bndIndices, reducibleBndData);
	auto reducibleBands = std::make_shared<elephon::ElectronicStructure::ElectronicBands>();
	auto gridNoSym = std::make_shared<elephon::LatticeStructure::RegularSymmetricGrid>();
	elephon::LatticeStructure::Symmetry idSym;
	idSym.set_reciprocal_space_sym();
	gridNoSym->initialize( 	bands->get_grid().get_grid_dim(),
							bands->get_grid().get_grid_prec(),
							bands->get_grid().get_grid_shift(),
							idSym,
							bands->get_grid().get_lattice()	);
	reducibleBands->initialize(bndIndices.size(), 0, reducibleBndData, *gridNoSym);

	auto tetra = std::make_shared<elephon::LatticeStructure::TetrahedraGrid>();
	tetra->initialize(gridNoSym);

	auto tetraSurf = std::make_shared< elephon::ElectronicStructure::TetrahedraIsosurface>();
	tetraSurf->initialize(tetra, reducibleBands, {0.0});

	tetraSurf->write_polydata_file( surfDataFile.string(), 0);
	boost::filesystem::remove(surfDataFile);
//	plot_iso_surface(tetraSurf, 0);


	auto tetraIrred = std::make_shared<elephon::LatticeStructure::TetrahedraGrid>();
	tetraIrred->initialize(std::make_shared<elephon::LatticeStructure::RegularSymmetricGrid>(bands->get_grid()));
	auto tetraSurfIrred = std::make_shared< elephon::ElectronicStructure::TetrahedraIsosurface>();
	tetraSurfIrred->initialize(tetraIrred, bands, {0.0});

	std::vector<double> dosRef;
	bands->compute_DOS_tetra( tetraIrred, {0}, dosRef);

	double dosReducible = 0;
	double dosIrreducible = 0;
	std::vector<double> kIso;
	std::vector<double> kWeights;
	for (int ibnd = 0 ; ibnd < bands->get_nBnd() ; ++ibnd)
	{
		tetraSurf->get_irreducible_iso_vector_integration_weights(0, ibnd, kIso, kWeights);
		for ( auto w : kWeights)
			dosReducible += w;
		tetraSurfIrred->get_irreducible_iso_vector_integration_weights(0, ibnd, kIso, kWeights);
		for ( auto w : kWeights)
			dosIrreducible += w;
	}
	BOOST_CHECK_CLOSE(dosRef[0], dosIrreducible, 1e-5);
	BOOST_CHECK_CLOSE(dosRef[0], dosReducible, 1e-5);
}

BOOST_AUTO_TEST_CASE( TetrahedraIsosurface_Al_sc444 )
{
	elephon::test::fixtures::MockStartup ms;
	auto rootDir = ms.get_data_for_testing_dir() / "Al" / "vasp" / "fcc_primitive" / "sc_4x4x4";
	auto surfDataFile = rootDir / "surface.vtp";
	boost::filesystem::remove(surfDataFile);
	elephon::test::fixtures::DataLoader dl;
	auto resHndler = dl.create_resource_handler(std::string()+
			"root_dir = "+rootDir.string()
			);
	auto tetraSurf = resHndler->get_tetrahedra_isosurface_fft();

	tetraSurf->write_polydata_file( surfDataFile.string(), 0);
	boost::filesystem::remove(surfDataFile);
//	plot_iso_surface(tetraSurf, 0);
}
