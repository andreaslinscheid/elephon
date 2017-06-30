set(elephon_LD_LIBS "")
set(elephon_INCLUDE "")
set(elephon_DEFINE "")
set(elephon_OPTS "")

# The projects own libraries
# %%%% HACK (!!!) %%%%
# ${subdirlist} appears twice to resolve cross dependencies of symbols
# the linker should stop searching once it hits the first, so this should work without causing harm ..
list(APPEND elephon_LD_LIBS  ${subdirlist} ${subdirlist} ${subdirlist})
list(APPEND elephon_INCLUDE ${elephon_src_tree})

# BOOST
find_package(Boost 1.36.0 COMPONENTS program_options unit_test_framework regex filesystem system REQUIRED)
list(APPEND elephon_LD_LIBS ${Boost_LIBRARIES} )
list(APPEND elephon_INCLUDE ${Boost_INCLUDE_DIRS})

# VTK 
find_package(VTK 6.0 REQUIRED NO_MODULE)
include( ${VTK_USE_FILE} )
list(APPEND elephon_LD_LIBS ${VTK_LIBRARIES} )

# FFTW3
list(APPEND elephon_LD_LIBS  "-lfftw3")

# LAPACKE
find_package( LAPACKE REQUIRED )
list(APPEND elephon_LD_LIBS ${LAPACKE_LIBRARIES} )

# BLAS
if(NOT APPLE)
  set(BLAS "Atlas" CACHE STRING "Selected BLAS library")
  set_property(CACHE BLAS PROPERTY STRINGS "Atlas;Open;MKL")
  if(BLAS STREQUAL "Atlas" OR BLAS STREQUAL "atlas")
    find_package(Atlas REQUIRED)
    list(APPEND elephon_INCLUDE PUBLIC ${Atlas_INCLUDE_DIR})
    list(APPEND elephon_LD_LIBS PUBLIC ${Atlas_LIBRARIES})
  elseif(BLAS STREQUAL "Open" OR BLAS STREQUAL "open")
    find_package(OpenBLAS REQUIRED)
    list(APPEND elephon_INCLUDE PUBLIC ${OpenBLAS_INCLUDE_DIR})
    list(APPEND elephon_LD_LIBS PUBLIC ${OpenBLAS_LIB})
  elseif(BLAS STREQUAL "MKL" OR BLAS STREQUAL "mkl")
    find_package(MKL REQUIRED)
    list(APPEND elephon_INCLUDE PUBLIC ${MKL_INCLUDE_DIR})
    list(APPEND elephon_LD_LIBS PUBLIC ${MKL_LIBRARIES})
    list(APPEND elephon_OPTS PUBLIC -DUSE_MKL)
  endif()
elseif(APPLE)
  find_package(vecLib REQUIRED)
  list(APPEND elephon_INCLUDE PUBLIC ${vecLib_INCLUDE_DIR})
  list(APPEND elephon_LD_LIBS PUBLIC ${vecLib_LINKER_LIBS})

  if(VECLIB_FOUND)
    if(NOT vecLib_INCLUDE_DIR MATCHES "^/System/Library/Frameworks/vecLib.framework.*")
      list(APPEND elephon_OPTS PUBLIC -DUSE_ACCELERATE)
    endif()
  endif()
endif()

