cmake_minimum_required(VERSION 3.28)
project(SoftwareRenderer)

include(FetchContent)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/$<CONFIGURATION>")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/$<CONFIGURATION>")
set(SOURCE_FILE
        src/Resources/tgaimage.h
        src/Resources/tgaimage.cpp
        src/Resources/Mesh.cpp src/Resources/Mesh.h
        src/Utilities/Math/geometry.h
        src/Utilities/Math/geometry.cpp
        src/Renderer/SoftwareRenderer/my_gl.hpp
        src/Renderer/SoftwareRenderer/my_gl.cpp
        src/Utilities/Interfaces/IShader.hpp
        src/Utilities/CycleTimer.hpp
        src/main.cpp
        src/Utilities/Adapters/TGAtoSDLAdapter.hpp
        src/Renderer/SoftwareRenderer/render.cpp
        src/Renderer/SoftwareRenderer/render.hpp
        src/pch.h
        src/Renderer/SoftwareRenderer/shaders/PhongShadowShader.h
        src/Renderer/SoftwareRenderer/shaders/DepthShader.h
        src/Renderer/SoftwareRenderer/shaders/DepthShaderImage.h
        src/Renderer/SoftwareRenderer/shaders/GouraudShader.h
        src/Utilities/RectBuffer.h
        src/engine.cpp
        src/engine.h
        src/Utilities/Interfaces/IRenderer.h
        src/Utilities/Math/math.cpp
        src/Utilities/Math/math.h
        src/Resources/Material.cpp
        src/Resources/Material.h
        src/Resources/Scene.cpp
        src/Resources/Scene.h
        src/Resources/Camera.cpp
        src/Resources/Camera.h
        src/OS/Window.cpp
        src/OS/Window.h
        src/Resources/Model.cpp
        src/Resources/Model.h)
add_compile_options( )
add_executable(${PROJECT_NAME} ${SOURCE_FILE})

if(CMAKE_CXX_COMPILER_ID MATCHES Clang)
    add_compile_options(-mavx2)
    # For Address Sanitizer, only available on Clang
    target_compile_options(${PROJECT_NAME} PRIVATE
            $<$<CONFIG:Release>: -ffast-math;-O2;-mavx2;-Rpass=loop-vectorize;-Rpass-missed=loop-vectorize>
            $<$<CONFIG:Debug>: -Wall;-fsanitize=leak;-fno-omit-frame-pointer;-fsanitize-blacklist=blacklist.txt>)
elseif(CMAKE_CXX_COMPILER_ID MATCHES GNU)
    # Auto vectorization info
    target_compile_options(${PROJECT_NAME} PRIVATE
            $<$<CONFIG:Release>: -O2;-mavx2;-fopt-info-vec;-fopt-info-vec-missed;-fopenmp;-ffast-math>
            $<$<CONFIG:RelWithDebInfo>:  -O2;-mavx2;-fopt-info-vec;-fopt-info-vec-missed;-ffast-math>
            $<$<CONFIG:Debug>:-Wall;-fopenmp>)
elseif(CMAKE_CXX_COMPILER_ID MATCHES MSVC)
    # /analyze: static memory analysis
    # /fsanitize=address: microsoft ver. of AddressSanitizer
    # Requires AddressSanitizer, will crash if its not enabled.
    # /Qvec-report:2: Auto vectorization info
    target_compile_options(${PROJECT_NAME} PRIVATE
            $<$<CONFIG:Release>: /O2;/arch:AVX2;/Qvec-report:2>
            $<$<CONFIG:Debug>:/analyze;/openmp:llvm>
            $<$<CONFIG:RelWithDebInfo>:/arch:AVX2;/Qvec-report:2;/openmp:llvm>)
endif()

# Copy the directory from source code tree to executable tree
file(COPY ./obj DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY blacklist.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
target_precompile_headers(${PROJECT_NAME} PUBLIC src/pch.h)

# --- LIBRARIES ---
# Add OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    target_link_libraries(SoftwareRenderer PUBLIC OpenMP::OpenMP_CXX)
endif()

# Add SDL
FetchContent_Declare(
        SDL
        GIT_REPOSITORY https://github.com/libsdl-org/SDL
        GIT_TAG         a962f40 #3.4.0
)
FetchContent_MakeAvailable(SDL)
#add_subdirectory(third_party/SDL EXCLUDE_FROM_ALL)
#add_compile_options(-DSDL_X11_XCURSOR=OFF)
target_link_libraries(SoftwareRenderer PUBLIC SDL3::SDL3)

# Add imgui
# imgui is only available as source code and does not have their own CMake config
# 1. Setup paths
set(IMGUI_DIR "${CMAKE_SOURCE_DIR}/src/third_party/imgui")

# 3. Add ImGui sources DIRECTLY to your executable
add_library(imgui STATIC
        "${IMGUI_DIR}/imgui.cpp"
#        "${IMGUI_DIR}/imgui_demo.cpp"
        "${IMGUI_DIR}/imgui_draw.cpp"
        "${IMGUI_DIR}/imgui_tables.cpp"
        "${IMGUI_DIR}/imgui_widgets.cpp"
        "${IMGUI_DIR}/backends/imgui_impl_sdlrenderer3.cpp"
        "${IMGUI_DIR}/backends/imgui_impl_sdl3.cpp"
)

# 4. Add the include directory so #include "imgui.h" works
target_include_directories(imgui PUBLIC "${IMGUI_DIR}" "${IMGUI_DIR}/backends")
target_link_libraries(imgui PUBLIC SDL3::SDL3) #imgui needs SDL3 as well since we are using their implementation functions
target_link_libraries(${PROJECT_NAME} PRIVATE imgui)