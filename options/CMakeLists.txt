cmake_minimum_required(VERSION 3.16)

project(options)
include_directories(include)

add_executable(options_explicit src/test_explicit.cpp)
add_executable(options_implicit src/test_implicit.cpp)

find_package(GSL REQUIRED)
target_link_libraries(options_explicit GSL::gsl GSL::gslcblas)
target_link_libraries(options_implicit GSL::gsl GSL::gslcblas)
