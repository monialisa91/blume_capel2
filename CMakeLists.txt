cmake_minimum_required(VERSION 3.15)
project(blume_capel)

set(CMAKE_CXX_STANDARD 14)

add_executable(blume_capel main.cpp metropolis.h wolff.h observables.h)
target_link_libraries(blume_capel armadillo)