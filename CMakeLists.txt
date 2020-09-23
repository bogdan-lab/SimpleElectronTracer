cmake_minimum_required(VERSION 3.0)

PROJECT(ElectronTracer)

#include_directories( /home/bogdan/gmsh_sdk/gmsh/api/)

set (CMAKE_VERBOSE_MAKEFILE ON)

set (CMAKE_CXX_STANDARD 17)

add_executable( electron_tracer
    src/electron_tracer.cpp
    src/reflector.cpp
    src/particle.cpp
    src/surface.cpp
    src/utils.cpp
)

#target_link_libraries(gmsh_reader /home/bogdan/gmsh_sdk/lib64/libgmsh.so)
target_link_libraries(electron_tracer pthread)
