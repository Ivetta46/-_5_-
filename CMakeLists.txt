cmake_minimum_required(VERSION 3.20)
project(__1)

set(CMAKE_CXX_STANDARD 11)

add_executable(__1 solve.cpp localisation.cpp polynomial.cpp polynomial.h)
