cmake_minimum_required(VERSION 2.8)

project(elephon)

if (CMAKE_VERSION VERSION_LESS "3.1")
  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set (CMAKE_CXX_FLAGS "--std=gnu++11 ${CMAKE_CXX_FLAGS}")
  endif ()
else ()
  set (CMAKE_CXX_STANDARD 11)
endif ()

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

find_package(Boost 1.36.0 COMPONENTS program_options unit_test_framework REQUIRED)
find_package(VTK 6.0 REQUIRED NO_MODULE)
include(${VTK_USE_FILE})


target_link_libraries(elephon PUBLIC
  ${Boost_LIBRARIES} ${subdirlist} ${VTK_LIBRARIES} 
)

include(CTest)

foreach(testSrc ${sources_test})
        #Extract the filename without an extension (NAME_WE)
        get_filename_component(testName ${testSrc} NAME_WE)

        #Add compile target
        add_executable(${testName} ${testSrc})

        #link to Boost libraries AND your targets and dependencies
        target_link_libraries(${testName} ${Boost_LIBRARIES} ${subdirlist} ${VTK_LIBRARIES})

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
