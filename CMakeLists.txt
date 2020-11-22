cmake_minimum_required(VERSION 3.16)
project(DendriteV3)

#include_directories(/home/xidad/libs/eigen-3.3.7)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -fopenmp")

set(CMAKE_VERBOSE_MAKEFILE ON)

add_executable(DendriteV3 MyTypedefs.h main.cpp Domain.cpp Domain.h BoundingRect.cpp BoundingRect.h Solver.cpp Solver.h)