cmake_minimum_required(VERSION 3.1)

project(ElectronTracer
	LANGUAGES CXX)

option(BUILD_TESTS OFF)
option(BUILD_BENCHMARKS OFF)

include_directories(include/)

if(BUILD_BENCHMARKS)
    set(CMAKE_BUILD_TYPE Release)
endif()


add_subdirectory(src)
if(BUILD_TESTS)
        enable_testing()
        add_subdirectory(tests)
endif()

if(BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()

