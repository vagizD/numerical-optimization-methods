cmake_minimum_required(VERSION 3.21)
project(grobner-bases)

set(CMAKE_CXX_STANDARD 23)

include_directories(include)
include_directories(src)

set(SOURCE_FILES
        include/algebra.hpp
        include/backtracking.hpp
        include/constants.hpp
        include/monomial.hpp
        include/ordering.hpp
        include/polynomial.hpp
        include/sturm.hpp
        include/system.hpp
        include/sturm.hpp src/main.cpp
)

add_executable(grobner-bases ${SOURCE_FILES})
