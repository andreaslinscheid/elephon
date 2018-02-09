
#define BOOST_TEST_MODULE "Unit Tests for Elephon"
#include <boost/test/included/unit_test.hpp>

// include all test suits for namespace Algorithms
#include "Algorithms/test_FFTInterface.cpp"
#include "Algorithms/test_LinearAlgebraInterface.cpp"
#include "Algorithms/test_TrilinearInterpolation.cpp"
#include "Algorithms/test_CubeSplineInterpolation.cpp"
#include "Algorithms/test_rotation_matrix_connecting_unit_vectors.cpp"
#include "Algorithms/test_SphereIntegrator.cpp"

// include all test suits for namespace AtomicSite
#include "AtomicSite/test_SphericalHarmonicExpansion.cpp"
#include "AtomicSite/test_WignerDMatrix.cpp"
#include "AtomicSite/test_EulerAngles.cpp"

// include all test suits for namespace Auxillary
#include "Auxillary/test_CubicPolynomial.cpp"

// include all test suits for namespace ElectronicStructure
#include "ElectronicStructure/test_BandStructureAnalysis.cpp"
#include "ElectronicStructure/test_ElectronicBands.cpp"
#include "ElectronicStructure/test_FermiSurface.cpp"
#include "ElectronicStructure/test_GradientFFTReciprocalGrid.cpp"
#include "ElectronicStructure/test_LocalDensityOfStates.cpp"
#include "ElectronicStructure/test_NestingFunction.cpp"
#include "ElectronicStructure/test_TetrahedraIsosurface.cpp"
#include "ElectronicStructure/test_Wavefunctions.cpp"

// include all test suits for namespace IOMethods
#include "IOMethods/test_BuildFolderStructure.cpp"
#include "IOMethods/test_Input.cpp"
#include "IOMethods/test_KPath.cpp"
#include "IOMethods/test_ReadVASPElephonData.cpp"
#include "IOMethods/test_ReadVASPLocpot.cpp"
#include "IOMethods/test_ReadVASPPoscar.cpp"
#include "IOMethods/test_ReadVASPSymmetries.cpp"
#include "IOMethods/test_ReadVASPWaveFunction.cpp"
#include "IOMethods/test_ReadVASPxmlFile.cpp"
#include "IOMethods/test_WriteVASPRealSpaceData.cpp"

// include all test suits for namespace LatticeStructure
#include "LatticeStructure/test_DataRegularGrid.cpp"
#include "LatticeStructure/test_LatticeModule.cpp"
#include "LatticeStructure/test_RegularBareGrid.cpp"
#include "LatticeStructure/test_RegularSymmetricGrid.cpp"
#include "LatticeStructure/test_Symmetry.cpp"
#include "LatticeStructure/test_SymmetryReduction.cpp"
#include "LatticeStructure/test_TetrahedraGrid.cpp"
#include "LatticeStructure/test_UnitCell.cpp"
#include "LatticeStructure/test_PrimitiveToSupercellConnection.cpp"

// include all test suits for namespace PhononStructure
#include "PhononStructure/test_AlphaSquaredF.cpp"
#include "PhononStructure/test_DisplacementPotential.cpp"
#include "PhononStructure/test_ElectronPhononCoupling.cpp"
#include "PhononStructure/test_ForceConstantMatrix.cpp"
#include "PhononStructure/test_Phonon.cpp"
#include "PhononStructure/test_PhononGrid.cpp"
#include "PhononStructure/test_PotentialChangeIrredDisplacement.cpp"

// include all test suits for namespace symmetry
#include "symmetry/test_SymmetryOperation.cpp"
