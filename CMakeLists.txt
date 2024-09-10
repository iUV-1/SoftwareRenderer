cmake_minimum_required(VERSION 3.28)
project(SoftwareRenderer)

set(CMAKE_CXX_STANDARD 20)

add_executable(SoftwareRenderer main.cpp
        tgaimage.h
        tgaimage.cpp
        model.cpp
        model.h
        geometry.h)
