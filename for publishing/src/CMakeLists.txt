cmake_minimum_required(VERSION 3.28)
project(SoftwareRenderer)

set(CMAKE_CXX_STANDARD 20)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libstdc++ -static-libgcc")

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
