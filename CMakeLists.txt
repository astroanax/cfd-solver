cmake_minimum_required(VERSION 3.10)
project(euler)

# ensure that the C++ 2017 standard is used
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(euler euler.cpp)