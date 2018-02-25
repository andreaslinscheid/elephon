cmake_minimum_required(VERSION 3.3)

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

list(REMOVE_ITEM VTK_LIBRARIES "vtkproj4")

#Break up the test sources from test fixtures
set(sources_test)
file(GLOB_RECURSE testModules ../test/*.cpp)
foreach(testSrc ${testModules})
	if ( NOT (${testSrc} MATCHES "test_main.cpp") )
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
list(INSERT elephon_LD_LIBS 0 ${elephonLibraryNames})
list(INSERT elephon_INCLUDE 0 ${elephon_src_tree})

#add a library with the tests
add_library( libTest ${sources_test} )

target_compile_options(libTest PUBLIC ${elephon_OPTS} )
target_include_directories(libTest PUBLIC ${elephon_src_tree} "${PROJECT_SOURCE_DIR}/../test/" ${elephon_INCLUDE})

add_executable(elephon main.cpp)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -ffast-math -mtune=generic -DBOOST_DISABLE_ASSERTS")


target_compile_options(elephon PUBLIC ${elephon_OPTS} )
target_include_directories(elephon PUBLIC ${elephon_INCLUDE})

target_link_libraries(elephon PUBLIC -Wl,--start-group ${elephon_LD_LIBS} -Wl,--end-group )

include(CTest)

#Add compile target
add_executable(test_main "${PROJECT_SOURCE_DIR}/../test/test_main.cpp")

#link to Boost libraries AND your targets and dependencies
target_link_libraries(test_main PUBLIC libTest -Wl,--start-group ${elephon_LD_LIBS} -Wl,--end-group )

#I like to move testing binaries into a test/bin directory
set_target_properties(test_main PROPERTIES 
    RUNTIME_OUTPUT_DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR}/../test/bin)
target_compile_options(test_main PUBLIC ${elephon_OPTS})
target_include_directories(test_main PUBLIC ${elephon_INCLUDE} "${PROJECT_SOURCE_DIR}/../test/")

SET(BOOST_TEST_LOG_LEVEL=test_suite)

#Finally add it to test execution - 
#Notice the WORKING_DIRECTORY and COMMAND
add_test(NAME test_main 
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../test/bin 
         COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../test/bin/test_main
         ARGUMENT --log_level=test_suite)
