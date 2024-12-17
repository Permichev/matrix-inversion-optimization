cmake_minimum_required(VERSION 3.15)
project(lab2)

set(CMAKE_CXX_STANDARD 17)

# Флаги оптимизации (при желании можно добавить или убрать некоторые)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Ofast -march=native -ffast-math")


add_executable(lab2 invert_matrix.cpp)
