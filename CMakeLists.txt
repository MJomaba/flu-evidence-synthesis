project (flu-evidence-synthesis)

cmake_minimum_required(VERSION 2.8)

SET (CMAKE_CXX_FLAGS "-O3 -std=c++11 -Wall")
#SET (CMAKE_CXX_FLAGS "-O2 -pipe -ffast-math -fomit-frame-pointer -Wfatal-errors -std=gnu++0x")

# Can be used by youcompleteme to lookup compile information
SET (CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(GSL REQUIRED)

set(BOOST_LIBS program_options) 
find_package(Boost COMPONENTS ${BOOST_LIBS} REQUIRED)


SET (CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)
add_executable (flu-evidence-synthesis src/main.cc)
target_link_libraries(flu-evidence-synthesis ${GSL_LIBRARIES};${Boost_LIBRARIES})
