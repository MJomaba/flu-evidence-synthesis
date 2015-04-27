project (flu-evidence-synthesis)

cmake_minimum_required(VERSION 2.8)

find_package(GSL REQUIRED)

SET (CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)
add_executable (flu-evidence-synthesis SRC/main.c)
target_link_libraries(flu-evidence-synthesis ${GSL_LIBRARIES})
