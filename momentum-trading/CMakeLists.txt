cmake_minimum_required(VERSION 3.28)

project(momentum_trading)
include_directories(include)

set(CMAKE_CXX_STANDARD 17)

add_executable(momentum_trading
        src/model.cpp)
