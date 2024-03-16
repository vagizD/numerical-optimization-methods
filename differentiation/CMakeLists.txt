cmake_minimum_required(VERSION 3.21)
project(differentiation)

set(CMAKE_CXX_STANDARD 20)

set(SOURCE_FILES
        src/aad.cpp
        src/tests.cpp
)

include_directories(include)
add_executable(diff ${SOURCE_FILES})
