set(elephon_LD_LIBS "")
set(elephon_INCLUDE "")
set(elephon_DEFINE "")
set(elephon_OPTS "")

include(../cmake/Macros.cmake)

# BOOST
prefer_static() 
find_package(Boost 1.56 COMPONENTS unit_test_framework regex filesystem system REQUIRED)
list(APPEND elephon_LD_LIBS ${Boost_LIBRARIES} )
list(APPEND elephon_INCLUDE ${Boost_INCLUDE_DIRS})

# VTK 
prefer_static() 
find_package(VTK REQUIRED NO_MODULE)
include( ${VTK_USE_FILE} )
list(APPEND elephon_LD_LIBS ${VTK_LIBRARIES} )

# FFTW3
prefer_static() 
find_package( FFTW3 )
list(APPEND elephon_INCLUDE ${FFTW3_INCLUDES})
list(APPEND elephon_LD_LIBS ${FFTW3_LIBRARIES})

# BLAS
prefer_static() 
if(NOT APPLE)
  set(BLAS "Open" CACHE STRING "Selected BLAS library")
  set_property(CACHE BLAS PROPERTY STRINGS "Atlas;Open;MKL")
  if(BLAS STREQUAL "Atlas" OR BLAS STREQUAL "atlas")
    find_package(Atlas REQUIRED)
    list(APPEND elephon_INCLUDE ${Atlas_INCLUDE_DIR})
    list(APPEND elephon_LD_LIBS ${Atlas_LIBRARIES})
    list(APPEND elephon_OPTS -DUSE_ATLAS)
  elseif(BLAS STREQUAL "Open" OR BLAS STREQUAL "open")
    find_package(OpenBLAS REQUIRED)
    list(APPEND elephon_INCLUDE ${OpenBLAS_INCLUDE_DIR})
    list(APPEND elephon_LD_LIBS ${OpenBLAS_LIB})
  elseif(BLAS STREQUAL "MKL" OR BLAS STREQUAL "mkl")
    find_package(MKL REQUIRED)
    list(APPEND elephon_INCLUDE ${MKL_INCLUDE_DIR})
    list(APPEND elephon_LD_LIBS ${MKL_LIBRARIES})
    list(APPEND elephon_OPTS -DUSE_MKL)
  endif()
elseif(APPLE)
  find_package(vecLib REQUIRED)
  list(APPEND elephon_INCLUDE ${vecLib_INCLUDE_DIR})
  list(APPEND elephon_LD_LIBS ${vecLib_LINKER_LIBS})

  if(VECLIB_FOUND)
    if(NOT vecLib_INCLUDE_DIR MATCHES "^/System/Library/Frameworks/vecLib.framework.*")
      list(APPEND elephon_OPTS PUBLIC -DUSE_ACCELERATE)
    endif()
  endif()
endif()

# OPENMP
prefer_static() 
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
