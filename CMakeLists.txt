cmake_minimum_required(VERSION 3.16)
project(DendriteV3)

set(CMAKE_CUDA_STANDARD 14)

#include_directories(/home/xidad/libs/eigen-3.3.7)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_CUDA_ARCHITECTURES 50)

add_executable(DendriteV3 MyTypedefs.h main.cpp Domain.cpp Domain.h BoundingRect.cpp BoundingRect.h Solver.cpp Solver.h)