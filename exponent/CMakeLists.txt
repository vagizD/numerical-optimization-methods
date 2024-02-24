cmake_minimum_required(VERSION 3.21)
project(exponent)

set(CMAKE_CXX_STANDARD 20)

find_package(GSL::gsl CONFIG REQUIRED)
add_compile_options(-lcblas -lm -lgsl)
add_link_options(-lcblas -lm -lgsl -L/usr/local/lib)

set(SOURCE_FILES
    test-exp.cpp
    constants.hpp
    poly.hpp
    exp.hpp
    test.hpp
)

add_executable(exponent ${SOURCE_FILES})