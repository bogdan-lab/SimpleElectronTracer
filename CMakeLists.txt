cmake_minimum_required(VERSION 3.19)

project(ElectronTracer
	LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(BUILD_TESTS OFF)
option(BUILD_BENCHMARKS OFF)
option(BUILD_SANITIZE OFF)

include_directories(include/)

#WARNINGS
add_library(full_set_warnings INTERFACE)
set(WARNINGS_AS_ERRORS FALSE)
include(./cmake/CompilerWarnings.cmake)
set_project_warnings(full_set_warnings)

#CONAN
find_program(CONAN conan)
if (NOT CONAN)
    message(FATAL_ERROR "Needs conan package manager: please visit https://conan.io and install")
endif ()
include(./cmake/conan.cmake)
conan_check(VERSION 1.14.4 REQUIRED)
conan_cmake_run(REQUIRES
                nlohmann_json/3.9.1
                lyra/1.5.1
                fmt/7.1.0
        OPTIONS BASIC_SETUP CMAKE_TARGETS
        BUILD missing)


if(BUILD_BENCHMARKS)
    set(CMAKE_BUILD_TYPE Release)
endif()

if(BUILD_SANITIZE)
    set (CMAKE_BUILD_TYPE Debug)
    set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
    set (CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
endif()

add_subdirectory(src)
if(BUILD_TESTS)
        enable_testing()
        add_subdirectory(tests)
endif()

if(BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()

