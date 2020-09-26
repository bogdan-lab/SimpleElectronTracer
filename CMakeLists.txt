cmake_minimum_required(VERSION 3.1)

project(ElectronTracer
	LANGUAGES CXX)

set (CMAKE_VERBOSE_MAKEFILE ON)


add_executable( electron_tracer
    src/main.cpp
    src/reflector.cpp
    src/particle.cpp
    src/surface.cpp
    src/utils.cpp
)

