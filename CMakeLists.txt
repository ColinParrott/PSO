cmake_minimum_required(VERSION 3.12)
project(PSO)

set(CMAKE_CXX_STANDARD 14)

add_library(PSO SHARED src/PSO.cpp include/PSO.h)