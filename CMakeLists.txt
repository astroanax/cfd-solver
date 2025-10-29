cmake_minimum_required(VERSION 3.10)
project(euler)

# ensure that the C++ 2017 standard is used
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 1D Euler solver
add_executable(euler euler.cpp)

# 2D Euler solver
add_executable(euler2d euler2d.cpp)