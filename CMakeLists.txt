cmake_minimum_required(VERSION 3.28)
project(SoftwareRenderer)

set(CMAKE_CXX_STANDARD 20)
set(SOURCE_FILE main.cpp
        tgaimage.h
        tgaimage.cpp
        model.cpp
        model.h
        geometry.h
        my_gl.cpp
        my_gl.hpp
        my_gl.cpp
        geometry.cpp
        shaders.hpp)

add_link_options(-g)
add_compile_options(-Wall -g)

# Copy the directory from source code tree to executable tree
file(COPY ./obj DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(${PROJECT_NAME} ${SOURCE_FILE})
