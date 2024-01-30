cmake_minimum_required(VERSION 3.21)
project(NumMethods)

set(CMAKE_CXX_STANDARD 20)

add_executable(NumMethods test-exp.cpp constants.hpp exp.hpp)
