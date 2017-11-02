cmake_minimum_required(VERSION 2.8)

project(elephon)

set(elephon_src_tree ${PROJECT_SOURCE_DIR})

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../cmake/Modules)

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

include(../cmake/dependencies.cmake)

#Break up the test sources from test fixtures
set(sources_test)
set(fixtures_test)
file(GLOB_RECURSE testModules ../test/*.cpp)
foreach(testSrc ${testModules})
	if ( ${testSrc} MATCHES "/fixtures/" )
		list(APPEND fixtures_test ${testSrc})
	else()
		list(APPEND sources_test ${testSrc})
	endif()
endforeach(testSrc)

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

list(APPEND elephon_OPTS "-Wall" "-Wno-sign-compare" "-Wno-deprecated-declarations")

#add a library each
set(elephonLibraryNames "")
foreach(libElephon ${subdirlist})
	file(GLOB_RECURSE sources ${libElephon}/*.cpp ${libElephon}/*.h)
	file(GLOB_RECURSE check_not_header_only ${libElephon}/*.cpp)
	if( check_not_header_only )
		add_library( ${libElephon} ${sources} )
		target_compile_options(${libElephon} PUBLIC ${elephon_OPTS} )
		target_include_directories(${libElephon} PUBLIC ${elephon_src_tree} ${elephon_INCLUDE})
		list(APPEND elephonLibraryNames ${libElephon})
	endif()
endforeach(libElephon)

# The projects own libraries
# %%%% HACK (!!!) %%%%
# ${elephonLibraryNames} appears multiple times to resolve cross dependencies of symbols
# the linker should stop searching once it hits the first, so this should work without causing harm ..
list(INSERT elephon_LD_LIBS 0 ${elephonLibraryNames} ${elephonLibraryNames} ${elephonLibraryNames})
list(INSERT elephon_INCLUDE 0 ${elephon_src_tree})

#add a library with the tests
add_library( libTest ${fixtures_test} )

target_compile_options(libTest PUBLIC ${elephon_OPTS} )
target_include_directories(libTest PUBLIC ${elephon_src_tree} "${PROJECT_SOURCE_DIR}/../test/" ${elephon_INCLUDE})

add_executable(elephon main.cpp)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -ffast-math -mtune=generic")


target_compile_options(elephon PUBLIC ${elephon_OPTS} )
target_include_directories(elephon PUBLIC ${elephon_INCLUDE})

target_link_libraries(elephon PUBLIC ${elephon_LD_LIBS} )

include(CTest)

foreach(testSrc ${sources_test})
       	#Extract the filename without an extension (NAME_WE)
       	get_filename_component(testName ${testSrc} NAME_WE)

       	#Add compile target
       	add_executable(${testName} ${testSrc})

       	#link to Boost libraries AND your targets and dependencies
       	target_link_libraries(${testName} PUBLIC libTest ${elephon_LD_LIBS} )

       	#I like to move testing binaries into a test/bin directory
       	set_target_properties(${testName} PROPERTIES 
       	    RUNTIME_OUTPUT_DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR}/../test/bin)
       	target_compile_options(${testName} PUBLIC ${elephon_OPTS})
       	target_include_directories(${testName} PUBLIC ${elephon_INCLUDE} "${PROJECT_SOURCE_DIR}/../test/")

       	#Finally add it to test execution - 
       	#Notice the WORKING_DIRECTORY and COMMAND
       	add_test(NAME ${testName} 
       	         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../test/bin 
       	         COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../test/bin/${testName} )
endforeach(testSrc)
