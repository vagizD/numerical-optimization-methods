cmake_minimum_required(VERSION 3.16)

project(options-2)
include_directories(include)

add_executable(options_2_test src/test.cpp src/system_params.cpp
        include/utils.h)
