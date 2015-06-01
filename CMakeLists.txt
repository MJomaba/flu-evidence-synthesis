project (flu-evidence-synthesis)

cmake_minimum_required(VERSION 2.8)

SET (CMAKE_CXX_FLAGS "-O3 -std=c++11 -Wall")
#SET (CMAKE_CXX_FLAGS "-O2 -pipe -ffast-math -fomit-frame-pointer -Wfatal-errors -std=gnu++0x")

# Can be used by youcompleteme to lookup compile information
SET (CMAKE_EXPORT_COMPILE_COMMANDS ON)

SET( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} . )
include_directories( "/usr/local/include" )

find_package(GSL REQUIRED)

set(BOOST_LIBS regex filesystem program_options system iostreams thread) 
find_package(Boost COMPONENTS ${BOOST_LIBS} REQUIRED)

find_library( MONGO mongoclient )
if(NOT MONGO)
    message( FATAL_ERROR "Could not find mongoclient library, which is needed for JSON output" )
endif()
include_directories( "/usr/local/include/mongo" )
include_directories( "${CMAKE_SOURCE_DIR}/include" )
include_directories( "${CMAKE_SOURCE_DIR}/src" )

set(SRCS src/model.cc;src/io.cc;src/state.cc;src/trace.cc;src/data.cc;src/contacts.cc;src/vaccine.cc;src/proposal.cc)

SET (CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)
add_executable (flu-evidence-synthesis src/main.cc;${SRCS})
target_link_libraries(flu-evidence-synthesis ${MONGO};${GSL_LIBRARIES};${Boost_LIBRARIES})

add_executable (inference src/inference.cc;${SRCS})
target_link_libraries(inference ${MONGO};${GSL_LIBRARIES};${Boost_LIBRARIES})

#Documentation
set(DOXYFILE_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/doc")
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
include(UseDoxygen)

# Tests
# Includes Catch in the project:
include_directories(tests src ${COMMON_INCLUDES})
enable_testing(true)  # Enables unit-testing.

add_executable (test_fes tests/test_fes.cc;${SRCS})
target_link_libraries(test_fes ${MONGO};${GSL_LIBRARIES};${Boost_LIBRARIES})
add_test( test_fes bin/test_fes )
