cmake_minimum_required(VERSION 3.1)

project(ElectronTracer
	LANGUAGES CXX)

set (CMAKE_VERBOSE_MAKEFILE ON)

enable_testing()

add_subdirectory(src)
add_subdirectory(tests)

