cmake_minimum_required(VERSION 3.16)
project(satellite-model)

include_directories(includes)

set(CMAKE_CXX_STANDARD 20)

add_executable(satellite-model
        src/test.cpp
        src/observer.cpp
)

