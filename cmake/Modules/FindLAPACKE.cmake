include(CheckIncludeFile)
include(CheckFunctionExists)

function(_find_library _names _check_function _result)
    find_library(_lib
        NAMES ${_names}
        PATHS /usr/local/lib /usr/lib /usr/local/lib64 /usr/lib64
        ENV LD_LIBRARY_PATH
        )
    set(CMAKE_REQUIRED_LIBRARIES ${_lib})
    check_function_exists(${_check_function} _library_works)
    if(${_library_works})
        set(${_result} ${_lib} PARENT_SCOPE)
    endif()
    unset(_lib CACHE)
endfunction()

function(_find_include_dir _names _hints _result)
    find_path(_include_dir
        NAMES ${_names}
        HINTS ${_hints}
    )
    set(_all_include_files_work TRUE)
    foreach(_name ${_names})
        check_include_file(${_include_dir}/${_name} _include_file_works)
        if(NOT _include_file_works)
            set(_all_include_files_work FALSE)
        endif()
    endforeach()
    if(${_all_include_files_work})
        set(${_result} ${_include_dir} PARENT_SCOPE)
    endif()
    unset(_include_dir CACHE)
endfunction()

set(LAPACKE_FOUND FALSE)
set(LAPACKE_LIBRARIES "NOTFOUND")
set(LAPACKE_INCLUDE_DIR "NOTFOUND")

_find_library(lapacke LAPACKE_dgesv LAPACKE_LIBRARIES)
_find_include_dir(lapacke.h /usr LAPACKE_INCLUDE_DIR)

if(NOT "${LAPACKE_LIBRARIES}" MATCHES "NOTFOUND")
    if(NOT "${LAPACKE_INCLUDE_DIR}" MATCHES "NOTFOUND")
        set(LAPACKE_FOUND TRUE)
    endif()
endif()
