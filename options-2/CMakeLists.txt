cmake_minimum_required(VERSION 3.16)

project(options-2)
include_directories(include)
include_directories(test)

add_executable(
        options_2_test
        src/test.cpp src/system_params.cpp
)

add_compile_options(-fsanitize=address)
add_link_options(-fsanitize=address)

find_package(GSL REQUIRED)
target_link_libraries(options_2_test GSL::gsl GSL::gslcblas)
