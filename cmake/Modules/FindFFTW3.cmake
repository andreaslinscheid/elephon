# - Find FFTW3
# Find the native FFTW includes and library
#
#  FFTW3_INCLUDES    - where to find fftw3.h
#  FFTW3_LIBRARIES   - List of libraries when using FFTW.
#  FFTW3_FOUND       - True if FFTW found.

set(FFTW3_INCLUDE_SEARCH_PATHS
  /usr/include/
  $ENV{HPC_FFTW_INC}
)

set(FFTW3_LIB_SEARCH_PATHS
  /usr/lib/
  $ENV{HPC_FFTW_LIB}
)

if (FFTW3_INCLUDES)
  set (FFTW3_FIND_QUIETLY TRUE)
endif (FFTW3_INCLUDES)

find_path (FFTW3_INCLUDES NAMES fftw3.h PATHS ${FFTW3_INCLUDE_SEARCH_PATHS} )

find_library (FFTW3_LIBRARIES NAMES fftw3 PATHS ${FFTW3_LIB_SEARCH_PATHS})

# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (FFTW3 DEFAULT_MSG FFTW3_LIBRARIES FFTW3_INCLUDES)

mark_as_advanced (FFTW3_LIBRARIES FFTW3_INCLUDES)
