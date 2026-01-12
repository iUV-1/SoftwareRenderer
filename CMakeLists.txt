cmake_minimum_required(VERSION 3.28)
project(SoftwareRenderer)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/$<CONFIGURATION>")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/$<CONFIGURATION>")
set(SOURCE_FILE
        tgaimage.h tgaimage.cpp
        model.cpp model.h
        geometry.h geometry.cpp
        my_gl.hpp my_gl.cpp
        shaders.hpp
        CycleTimer.hpp
        main.cpp
        TGAtoSDLAdapter.cpp TGAtoSDLAdapter.hpp
        render.cpp render.hpp
        pch.h

        shaders/PhongShadowShader.h
        shaders/DepthShader.h
        shaders/DepthShaderImage.h)

#add_link_options(-g)
add_compile_options(-Wall -g -O2 -ffast-math )
add_executable(${PROJECT_NAME} ${SOURCE_FILE})

#add_compile_options(-fopt-info-vec
#        -fopt-info-vec-missed)
# Copy the directory from source code tree to executable tree
file(COPY ./obj DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
# 1. Find the OpenMP package
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    target_link_libraries(SoftwareRenderer PUBLIC OpenMP::OpenMP_CXX)
endif()

# Add SDL
add_subdirectory(third_party/SDL EXCLUDE_FROM_ALL)


target_precompile_headers(${PROJECT_NAME} PUBLIC pch.h)

target_link_libraries(SoftwareRenderer PUBLIC SDL3::SDL3)
