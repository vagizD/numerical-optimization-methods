cmake_minimum_required(VERSION 3.16)

project(options)
include_directories(include)

add_executable(options
        src/test.cpp
        include/constants.h
        include/stepper.h
        include/time_marshalling.h
)
