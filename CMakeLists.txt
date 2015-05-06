project (flu-evidence-synthesis)

cmake_minimum_required(VERSION 2.8)

SET (CMAKE_CXX_FLAGS "-O3 -std=c++11 -Wall")
#SET (CMAKE_CXX_FLAGS "-O2 -pipe -ffast-math -fomit-frame-pointer -Wfatal-errors -std=gnu++0x")

# Can be used by youcompleteme to lookup compile information
SET (CMAKE_EXPORT_COMPILE_COMMANDS ON)

SET( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} . )
include_directories( "/usr/local/include" )

find_package(GSL REQUIRED)

set(BOOST_LIBS filesystem program_options system iostreams) 
find_package(Boost COMPONENTS ${BOOST_LIBS} REQUIRED)

SET (CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)
add_executable (flu-evidence-synthesis src/main.cc;src/model.cc)
target_link_libraries(flu-evidence-synthesis ${GSL_LIBRARIES};${Boost_LIBRARIES})

#Documentation
set(DOXYFILE_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/doc")
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
include(UseDoxygen)

# Tests
# Includes Catch in the project:
include_directories(tests src ${COMMON_INCLUDES})
enable_testing(true)  # Enables unit-testing.

add_executable (test_fes tests/test_fes.cc)
target_link_libraries(test_fes ${GSL_LIBRARIES};${Boost_LIBRARIES})
add_test( test_fes bin/test_fes )
