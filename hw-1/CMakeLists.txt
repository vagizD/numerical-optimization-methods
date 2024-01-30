cmake_minimum_required(VERSION 3.21)
project(exponent)

set(CMAKE_CXX_STANDARD 20)

set(SOURCE_FILES
    test-exp.cpp
    constants.hpp
    exp.hpp
)

add_executable(exponent ${SOURCE_FILES})
