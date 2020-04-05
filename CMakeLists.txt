cmake_minimum_required(VERSION 3.0)

PROJECT(ElectronTracer)

#include_directories( /home/bogdan/gmsh_sdk/gmsh/api/)

set (CMAKE_VERBOSE_MAKEFILE ON)

set (CMAKE_CXX_STANDARD 17)
#set (CMAKE_C_STANDARD   11)
#set(CMAKE_INTERPROCEDURAL_OPTIMIZATION FALSE)

#set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}     -g -Wall  -Wextra -pedantic -Wno-sign-compare" )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}   -g -Wall  -Wextra -pedantic -Wno-sign-compare" )

#OPTION(SANITAIZE   "USE GCC sanitaizers" OFF)
#OPTION(VTUNE_BUILD "BUILD FOR VTUNE PROFILING" OFF)


#if(  ( CMAKE_BUILD_TYPE STREQUAL "Release") OR (NOT CMAKE_BUILD_TYPE) )
#    message("optimized build")
#    #set(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS}   -O3")
#    set(CMAKE_CXX_FLAGS   "${CMAKE_C_FLAGS}   -O3")
#    if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
#        set(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS}     -fopenmp -march=native ")
#        set(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS}   -fopenmp -march=native ")
#    elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
#	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -march=native" )
#	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas -Wno-overloaded-virtual -Wno-undefined-var-template -Wno-unused-function -Wno-unused-variable -Wno-inconsistent-missing-override")
#	set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -fopenmp -march=native" )
#    elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
#        if(VTUNE_BUILD)
#            # For Vtune need to have dynamic loaded librariries
#            set(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS}     -qopenmp -shared-intel -debug full -qopenmp-link dynamic -xHOST -parallel-source-info=2")
#            set(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS}   -qopenmp -shared-intel -debug full -qopenmp-link dynamic -xHOST -parallel-source-info=2")
#        else(VTUNE_BUILD)
#            set(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS}     -qopenmp -static-intel -xHOST")
#            set(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS}   -qopenmp -static-intel -xHOST")
#        endif(VTUNE_BUILD)
#    endif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    
#else( )
#    message("debug build")

#    if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
#        set(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS}      -fstack-protector -O0 -g3 -ggdb3 -fopenmp -mavx2") #-fopenmp
#        set(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS}    -fstack-protector -O0 -g3 -ggdb3 -fopenmp -mavx2") #-fopenmp
#        if(SANITAIZE)
#            message("Enable sanitaizers")
#            set(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS}     -fsanitize=address -fsanitize=bounds -fsanitize=bounds-strict") #-fopenmp
#            set(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS}   -fsanitize=address -fsanitize=bounds -fsanitize=bounds-strict ") #-fopenmp
#        endif(SANITAIZE)
#    else(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
#        set(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS}     -qopenmp -static-intel -O0")
#        set(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS}   -qopenmp -static-intel -O0")
#    endif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
#endif( )


add_executable( electron_tracer
    src/electron_tracer.cpp
)

#target_link_libraries(gmsh_reader /home/bogdan/gmsh_sdk/lib64/libgmsh.so)

install(TARGETS electron_tracer  RUNTIME DESTINATION ${ElectronTracer_SOURCE_DIR}/bin)
