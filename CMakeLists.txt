cmake_minimum_required(VERSION 3.28)
project(SoftwareRenderer)

set(CMAKE_CXX_STANDARD 20)
set(SOURCE_FILE main.cpp
        tgaimage.h tgaimage.cpp
        model.cpp model.h
        geometry.h geometry.cpp
        my_gl.hpp my_gl.cpp
        shaders.hpp)

add_link_options(-g)
add_compile_options(-Wall -g -O2 -ffast-math -fopt-info-vec
        -fopt-info-vec-missed)

# Copy the directory from source code tree to executable tree
file(COPY ./obj DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
# 1. Find the OpenMP package

add_executable(${PROJECT_NAME} ${SOURCE_FILE})
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    target_link_libraries(SoftwareRenderer PUBLIC OpenMP::OpenMP_CXX)
endif()