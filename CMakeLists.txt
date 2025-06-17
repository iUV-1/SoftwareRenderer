cmake_minimum_required(VERSION 3.28)
project(SoftwareRenderer)

set(CMAKE_CXX_STANDARD 20)
add_link_options(-g)
add_compile_options(-Wall -Wextra -O2 -g)

add_executable(SoftwareRenderer main.cpp
        tgaimage.h
        tgaimage.cpp
        model.cpp
        model.h
        geometry.h
        my_gl.cpp
        my_gl.hpp
        my_gl.cpp
        geometry.cpp)
