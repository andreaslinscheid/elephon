
#the usual place for intel libraries
set(INTEL_ROOT "/opt/intel" CACHE PATH )

find_path(MKL_ROOT include/mkl.h PATHS $ENV{MKLROOT} ${INTEL_ROOT}/mkl )

find_path(MKL_INCLUDE_DIR mkl.h PATHS ${MKL_ROOT} PATH_SUFFIXES include )

#
if(CMAKE_SIZEOF_VOID_P EQUAL 4)
  set(__path_suffixes lib lib/ia32)
else()
  set(__path_suffixes lib lib/intel64)
endif()
