#include <vector>
#include <cmath>
#include <iostream>
#include <limits>
#include <chrono>
#include <ctime>
#include <sstream>

#include "model.h"
#include "geometry.h"
#include "my_gl.hpp"
#include "tgaimage.h"

// this is due to gross coding from my end
// TODO: implement proper coding standards and untangle this mess
// Vec3f rasterize(Vec3f v, Matrix4x4f &m_viewport, Matrix4x4f const &m_proj, Matrix4x4f const &m_modelview) {
//     Matrix<float> homogonized = m_viewport*m_proj*m_modelview*homogonize(v);
//     Vec3f result = dehomogonize(homogonized);
//     // Round the result to apply to screen
//     result.x = std::round(result.x);
//     result.y = std::round(result.y);
//     result.z = std::round(result.z);
//
//     return result;
// }



float c;

Model *model = nullptr;

const int width = 800;
const int height = 800;

Vec3f light_dir = Vec3f(0.0, 0.0, 1.0);

struct GouraudShader: IShader {
    Vec3f varying_intensity;
    Matrix<float> vertex(int iface, int nthvert) override{
        Vec3f v = model->vert(iface, nthvert);
        Vec3f n = model->normal(iface, nthvert);
        // Cap at 0
        varying_intensity[nthvert] = std::max(0.f, n*light_dir);
        return Viewport*Projection*ModelView*homogonize(v);
    }
    // bar is barycentric of that vertex
    bool fragment(Vec3f bar, TGAColor &color) override {
        float intensity = varying_intensity*bar;
        color = TGAColor(255, 255, 255, 255) * intensity;
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
    if(argc < 2) {
        model = new Model("obj/african_head.obj");
    } else {
        model = new Model(argv[1]);
    }

    TGAImage tex_file(1024,1024,TGAImage::RGB);
    if(argc < 3) {
        tex_file.read_tga_file("obj/UV Grid.tga");
    } else {
        tex_file.read_tga_file(argv[2]);
    }

    tex_file.flip_vertically();

    // Setup zbuffer
    auto *zbuffer = new float[width * height];
    std::fill(zbuffer, zbuffer + width*height, -std::numeric_limits<float>::max()); // set every value in zbuffer to -inf

    TGAImage frame(width, height, TGAImage::RGB);

    // camera setting
    Vec3f eye(3, 2, 3);
    Vec3f cam(0, 0, 0);
    Vec3f up(0, 1, 0);

    // Setup GL
    LookAt(eye, cam, up);
    Project(5);

    GouraudShader shader = GouraudShader();
    light_dir.normalize();
    for (int i=0; i<model->nfaces(); ++i) {
        std::vector<int> face = model->face(i);
        std::vector<int> face_tex = model->face_tex(i);
        Vec3f screen_coords[3];
        Vec3f world_coords[3];
        Vec2f texture_coords[3];

        for (int j=0; j<3; ++j) {
            Vec3f v = model->vert(face[j]);
            world2screen(v, width, height, 255.0f);
            screen_coords[j] = rasterize(&shader, i, j);
            texture_coords[j] = model->texcoord(face_tex[j]);
        }
        // calculate normal
        // ^ is an overloaded operator that performs cross product calculation
        // world_coords[2] - world_coords[0] and the other are 2 vectors pointing from point
        // world_coords[0].
        Vec3f n = (screen_coords[2]-screen_coords[0])^(screen_coords[1]-screen_coords[0]);
        n.normalize();
        // calculate light intensity by dot product between normal and light vector
        float view_dir_intensity = eye*n;
        // back face culling
        //triangle(screen_coords, frame, zbuffer, TGAColor(intensity*255, intensity*255, intensity*255, 255));
        if (view_dir_intensity<1) {
            triangle(screen_coords, frame, zbuffer,
            tex_file, texture_coords, width, shader);
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

    delete[] zbuffer;
    delete model;
    return 0;
}