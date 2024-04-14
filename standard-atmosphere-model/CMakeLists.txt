cmake_minimum_required(VERSION 3.21)
project(standard-atmosphere-model)

set(CMAKE_CXX_STANDARD 20)

set(SOURCE_FILES
        src/sam.cpp
        src/tests.cpp
        src/observer.cpp
)

include_directories(include)
add_executable(sam ${SOURCE_FILES})
