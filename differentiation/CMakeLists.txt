cmake_minimum_required(VERSION 3.21)
project(differentiation)

set(CMAKE_CXX_STANDARD 20)

set(SOURCE_FILES
        src/diff.cpp
        src/stencil.cpp
        src/test_diff.cpp
)

add_executable(diff ${SOURCE_FILES})
