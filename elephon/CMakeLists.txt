cmake_minimum_required(VERSION 2.8)

project(elephon)

#Switch on C++11
if (CMAKE_VERSION VERSION_LESS "3.1")
  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set (CMAKE_CXX_FLAGS "--std=gnu++11 ${CMAKE_CXX_FLAGS}")
  endif ()
else ()
  set (CMAKE_CXX_STANDARD 11)
endif ()

#Easier error parsing
IF(CMAKE_COMPILER_IS_GNUCC)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fmessage-length=0")
ENDIF(CMAKE_COMPILER_IS_GNUCC)
IF(CMAKE_COMPILER_IS_GNUCXX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmessage-length=0")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

file(GLOB_RECURSE sources_test ../test/*.cpp)

#build a list of directories in elephon - this are the libraries we build
MACRO(SUBDIRLIST result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
      LIST(APPEND dirlist ${child})
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()
SUBDIRLIST( subdirlist ${CMAKE_CURRENT_SOURCE_DIR} )

#add a library each
foreach(libElephon ${subdirlist})
	file(GLOB_RECURSE sources ${libElephon}/*.cpp ${libElephon}/*.h)
	add_library( ${libElephon} ${sources} )
	target_include_directories(${libElephon} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
endforeach(libElephon)

add_executable(elephon main.cpp)

target_include_directories(elephon PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -ffast-math -march=native")

find_package(Boost 1.36.0 COMPONENTS program_options unit_test_framework filesystem system REQUIRED)
find_package(VTK 6.0 REQUIRED NO_MODULE)
include(${VTK_USE_FILE})

#HACK (!!!) - ${subdirlist} appears twice to resolve cross dependencies of symbols
# the linker should stop searching once it hits the first, so this should work without causing harm ..
target_link_libraries(elephon PUBLIC
  ${subdirlist} ${subdirlist} ${Boost_LIBRARIES} ${VTK_LIBRARIES} "-lfftw3"
)

include(CTest)

foreach(testSrc ${sources_test})
        #Extract the filename without an extension (NAME_WE)
        get_filename_component(testName ${testSrc} NAME_WE)

        #Add compile target
        add_executable(${testName} ${testSrc})

        #link to Boost libraries AND your targets and dependencies
	#HACK see before - ${subdirlist} appears twice
        target_link_libraries(${testName} ${subdirlist} ${subdirlist} ${Boost_LIBRARIES} ${VTK_LIBRARIES} "-lfftw3" )

        #I like to move testing binaries into a test/bin directory
        set_target_properties(${testName} PROPERTIES 
            RUNTIME_OUTPUT_DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR}/../test/bin)
	target_compile_options(${testName} PUBLIC -std=c++11 -Wall)
	target_include_directories(${testName} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

        #Finally add it to test execution - 
        #Notice the WORKING_DIRECTORY and COMMAND
        add_test(NAME ${testName} 
                 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../test/bin 
                 COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../test/bin/${testName} )
endforeach(testSrc)
