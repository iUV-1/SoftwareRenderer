#include <vector>
#include <iostream>
#include <limits>
#include <chrono>
#include <ctime>
#include <sstream>

#include "model.h"
#include "geometry.h"
#include "my_gl.hpp"
#include "tgaimage.h"

float c;

Model *model = nullptr;
TGAImage tex_file(1024,1024,TGAImage::RGB);
TGAImage normal_file(1024, 1024, TGAImage::RGB);

const int width = 800;
const int height = 800;

Vec3f light_dir = Vec3f(0.0, 0.0, 1.0);

struct GouraudShader: IShader {
    Vec3f varying_intensity; // intensity of a vertex
    Matrix<float> varying_uv = Matrix<float>(2, 3); // 2x3 matrix containing uv coordinate of 3 vertex (a trig)
    Matrix4x4f uniform_M; // Projection*ModelView
    Matrix4x4f uniform_MIT; // same as above but invert_transpose()
    Matrix4x4f uniform_view;

    Matrix<float> vertex(int iface, int nthvert) override{
        Vec3f v = model->vert(iface, nthvert);
        Vec3f n = model->normal(iface, nthvert);
        // Set the column of varying_uv to texture position in Vec2f
        varying_uv[0][nthvert] = model->texcoord(iface, nthvert).x;
        varying_uv[1][nthvert] = model->texcoord(iface, nthvert).y;
        // Cap at 0
        varying_intensity[nthvert] = std::max(0.f, n*light_dir);
        return uniform_view*homogonize(v);
    }
    // bar is the barycentric of that vertex
    bool fragment(Vec3f bar, TGAColor &color) override {
        // Cap at 0
        float intensity = std::max(0.f, varying_intensity*bar);
        // Convert barycentric vector to a matrix
        // NOTE: Somehow making a new variable is faster than making it inline?
        Matrix bary(bar); // 1x3 row matrix that represent a vector
        // Matrix<float> uv = varying_uv*Matrix(bar) <- Slower!
        Matrix<float> uv = varying_uv*bary; // 1x2 Matrix (Basically a Vec2f)

        // normal vector
        // Transform the normal vector to the eye space
        Vec3f n = dehomogonize(uniform_MIT*homogonize(model->normal(uv[0][0], uv[1][0]))).normalize();
        // light vector
        // Same as above
        Vec3f l = dehomogonize(uniform_M *homogonize(light_dir)).normalize();
        // omfg...
        //float intensity = std::max(0.f, n*l);
        TGAColor texColor = tex_file.get(uv[0][0] * tex_file.get_width(), uv[1][0] * tex_file.get_height());
        color = texColor * intensity;
        return false;
    }
};

Vec3f rasterize(GouraudShader *shader, int iface, int nthvert) {
    // Apply vertex shader
    Matrix<float> homogonized = shader->vertex(iface, nthvert);
    Vec3f result = dehomogonize(homogonized);
    // Round the result to apply to screen
    result.x = std::round(result.x);
    result.y = std::round(result.y);
    result.z = std::round(result.z);

    return result;
}

int main(int argc, char** argv) {
    // Time the render
    auto before = std::chrono::system_clock::now();

    if(argc < 2)
        model = new Model("obj/african_head.obj");
    else
        model = new Model(argv[1]);


    //TGAImage tex_file(1024,1024,TGAImage::RGB);
    if(argc < 3)
        tex_file.read_tga_file("obj/UV Grid.tga");
    else
        tex_file.read_tga_file(argv[2]);


    tex_file.flip_vertically();

    // TGAImage normal_file(1024, 1024, TGAImage::RGB);
    // if(argc < 4) {
    //     // TODO: ADD PLACEHOLDER NORMAL MAP
    //     normal_file.read_tga_file("placeholder");
    // } else {
    //     normal_file.read_tga_file(argv[4]);
    // }

    // Setup zbuffer
    float zbuffer[width*height];
    std::fill(zbuffer, zbuffer + width*height, -std::numeric_limits<float>::max()); // set every value in zbuffer to -inf

    TGAImage frame(width, height, TGAImage::RGB);

    // camera setting
    Vec3f eye(3, 2, 3);
    Vec3f cam(0, 0, 0);
    Vec3f up(0, 1, 0);

    // Setup GL
    LookAt(eye, cam, up);
    Project(5);
    SetViewport(width, height, 255.0f);

    GouraudShader shader = GouraudShader();
    shader.uniform_M = Projection*ModelView;
    shader.uniform_MIT = shader.uniform_M;
    shader.uniform_MIT.inverseTranspose();
    shader.uniform_view = Viewport*Projection*ModelView;
    light_dir.normalize();
    auto render = std::chrono::system_clock::now();
    for (int i=0; i<model->nfaces(); ++i) {
        Vec3f screen_coords[3];

        for (int j=0; j<3; ++j)
            screen_coords[j] = rasterize(&shader, i, j);

        // calculate normal
        // ^ is an overloaded operator that performs cross product calculation
        // world_coords[2] - world_coords[0] and the other are 2 vectors pointing from point
        // world_coords[0].
        Vec3f n = (screen_coords[2]-screen_coords[0])^(screen_coords[1]-screen_coords[0]);
        n.normalize();
        // calculate eye intensity by dot product between normal and eye vector
        float view_dir_intensity = eye*n;
        // back face culling
        if (view_dir_intensity<1) {
            triangle(screen_coords, frame, zbuffer, width, shader);
        }
    }
    // set origin to the bottom left corner
    frame.flip_vertically();
    // Get timing of the render
    auto now = std::chrono::system_clock::now();
    auto finish_time = std::chrono::system_clock::to_time_t(now);
    std::tm local_time = *std::localtime(&finish_time);
    std::stringstream sstream;
    sstream << "../output/" << local_time.tm_mon + 1 << "-" << local_time.tm_mday << "_"
            << local_time.tm_hour << "-" << local_time.tm_min << ".tga";
    frame.write_tga_file(sstream.str().c_str());

    // How long the render takes
    std::chrono::duration<double> diff = now - before;
    auto ms = duration_cast<std::chrono::milliseconds>(diff);
    std::cout << "Elapsed time: " << ms.count() << " ms" << std::endl;
    diff = now-render;
    ms = duration_cast<std::chrono::milliseconds>(diff);
    std::cout << "Render time: " << ms.count() << " ms" << std::endl;

    delete model;
    return 0;
}