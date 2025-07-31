//
// Created by iuv on 7/30/25.
//

#ifndef SHADERS_HPP
#define SHADERS_HPP
//
// Created by iuv on 7/30/25.
//

#include <iostream>
#include <utility>

#include "model.h"
#include "shaders.hpp"
#include "geometry.h"
#include "my_gl.hpp"

Vec3f light = Vec3f(-1.0, 1.0, 1.0).normalize();
extern Model *model;

extern TGAImage tex_file;
extern TGAImage normal_file;
extern TGAImage specular_file;

extern bool use_specular; // Use specular map
extern bool use_normal; // Use normal map

extern const float depth;
extern const int width;

struct GouraudShaderReference: IShader {
    Vec3f varying_intensity; // intensity of a vertex

    Matrix<float> vertex(int iface, int nthvert) override{
        Vec3f v = model->vert(iface, nthvert);
        Vec3f n = model->normal(iface, nthvert);
        // Set the column of varying_uv to texture position in Vec2f
        varying_uv.set_col(nthvert, model->texcoord(iface, nthvert));
        // Cap at 0
        varying_intensity[nthvert] = std::max(0.f, n*light);
        return Viewport*uniform_M*homogonize(v, 1.);
    }
    // bar is the barycentric of that vertex
    bool fragment(Vec3f bar, TGAColor &color) override {
        // Cap at 0
        float intensity = std::max(0.f, varying_intensity*bar);
        // Convert barycentric vector to a matrix
        // NOTE: Somehow making a new variable is faster than making it inline?
        Matrix<float> bary = Matrix(bar); // 1x3 row matrix that represent a vector
        // Matrix<float> uv = varying_uv*Matrix(bar) <- Slower!
        Matrix<float> uv = varying_uv*bary; // 1x2 Matrix (Basically a Vec2f)
        TGAColor texColor = tex_file.get(uv[0][0] * tex_file.get_width(), uv[1][0] * tex_file.get_height());
        color = texColor * intensity;
        return false;
    }
};

// Gouraud Shader uses vertex data to calculate light value
// However, this shader uses interpolation to compute the normal vector per pixel
// This is also called smooth shading
struct GouraudShader: IShader {
    Vec3f varying_intensity; // intensity of a vertex
    Matrix<float> varying_uv = Matrix<float>(2, 3); // 2x3 matrix containing uv coordinate of 3 vertex (a trig)
    Matrix4x4f uniform_M; // Projection*ModelView
    Matrix4x4f uniform_MIT; // same as above but invert_transpose()

    Matrix<float> vertex(int iface, int nthvert) override{
        //light.normalize();
        Vec3f v = model->vert(iface, nthvert);
        Vec3f n = model->normal(iface, nthvert);
        // Set the column of varying_uv to texture position in Vec2f
        varying_uv[0][nthvert] = model->texcoord(iface, nthvert).x;
        varying_uv[1][nthvert] = model->texcoord(iface, nthvert).y;
        // Cap at 0
        varying_intensity[nthvert] = std::max(0.f, n*light);
        return Viewport*Projection*ModelView*homogonize(v, 1.);
    }
    // bar is the barycentric of that vertex
    bool fragment(Vec3f bar, TGAColor &color) override {
        // Cap at 0
        // Convert barycentric vector to a matrix
        // NOTE: Somehow making a new variable is faster than making it inline?
        Matrix<float> bary = Matrix(bar); // 1x3 row matrix that represent a vector
        // Matrix<float> uv = varying_uv*Matrix(bar) <- Slower!
        Matrix<float> uv = varying_uv*bary; // 1x2 Matrix (Basically a Vec2f)

        Vec3f norm;
        if(false)
            norm = model->normal(uv[0][0], uv[1][0]);
        else {
            TGAColor normal_color = normal_file.get(uv[0][0] * normal_file.get_width(), uv[1][0] * normal_file.get_height());
            norm = Vec3f(normal_color.r, normal_color.g, normal_color.b);
        }
        // Transforming
        Vec3f l = dehomogonize(uniform_M  *homogonize(light, 0.f)).normalize();
        //l.z = -l.z;
        l.x = -l.x;
        l.y = -l.y;
        l.z = -l.z;
        Vec3f n = dehomogonize(uniform_MIT*homogonize(norm, 0.f)).normalize();
/*        if(varying_intensity*bar != n*l) {
            std::cout << "mismatched value" << std::endl;
        }*/
        float diff = std::max(0.f, n * l); // diffuse intensity value

        TGAColor texColor = tex_file.get(uv[0][0] * tex_file.get_width(), uv[1][0] * tex_file.get_height());
        color = texColor * diff;
        return false;
    }
};

// Phong Shader includes Phong reflection model and specular mapping
struct PhongShader: IShader {
    Matrix<float> varying_uv = Matrix<float>(2, 3); // 2x3 matrix containing uv coordinate of 3 vertex (a trig)
    Matrix4x4f uniform_M; // Projection*ModelView
    Matrix4x4f uniform_MIT; // same as above but invert_transpose()

    Matrix<float> vertex(int iface, int nthvert) override{
        Vec3f v = model->vert(iface, nthvert);
        // Set the column of varying_uv to texture position in Vec2f
        varying_uv.set_col(nthvert, model->texcoord(iface, nthvert));

        return Viewport*uniform_M*homogonize(v, 1.);
    }
    // bar is the barycentric of that vertex
    bool fragment(Vec3f bar, TGAColor &color) override {
        // Convert barycentric vector to a matrix
        // NOTE: Somehow making a new variable is faster than making it inline?
        Matrix<float> bary = Matrix(bar); // 1x3 row matrix that represent a vector
        // Matrix<float> uv = varying_uv*Matrix(bar) <- Slower!
        Matrix<float> mat_uv = varying_uv*bary; // 1x2 Matrix (Basically a Vec2f)
        Vec2f uv = Vec2f(mat_uv[0][0], mat_uv[1][0]);
        // Get the normal vector of that mesh based on the setting
        Vec3f norm;
        if(!use_normal)
            norm = model->normal(uv.u, uv.v);
        else {
            TGAColor normal_color = normal_file.get(uv.u * normal_file.get_width(), uv.v * normal_file.get_height());
            norm = Vec3f(normal_color.r, normal_color.g, normal_color.b);
        }
        /// Insanely costly calculations
        // Transform the normal vector to the eye space
        Vec3f n = dehomogonize(uniform_MIT*homogonize(norm, 0.f)).normalize();
        // Same as above
        Vec3f l = dehomogonize(uniform_M  *homogonize(light, 0.f)).normalize();
        l.z = -l.z; // I think this formula is meant to work for a different axis?
        l.y = -l.y;
        l.x = -l.x;
        Vec3f r = (n*(n*l*2.f) - l).normalize(); // reflection vector
        float diff = std::max(0.f, n * l); // diffuse intensity value
        // Specular
        float spec = 0.f;
        if(use_specular) {
            float spec_map_val = specular_file.get(uv.u * specular_file.get_width(), uv.v * specular_file.get_height()).r;
            spec = pow(std::max(r.z, 0.f), spec_map_val);
        } else {
            spec = std::max(r.z, 0.f);
        }
        TGAColor texColor = tex_file.get(uv.u * tex_file.get_width(), uv.v * tex_file.get_height());
        for (int i = 0; i < 3; i++) color[i] = std::min<float>(5 + texColor[i]*(diff + .6*spec), 255.f);
        color = texColor * diff;

        return false;
    }
};

// Like above but includes a shadow pass
struct PhongShaderShadow: IShader {
    Matrix3x3<float> varying_tri; // 3x3 matrix containing verticies of a trig
    Matrix4x4f uniform_Mshadow; // Shadow transformation
    float* depth_buffer;

    PhongShaderShadow(Matrix4x4f uniform_shadow, float* depth_buffer) : uniform_Mshadow(std::move(uniform_shadow)), depth_buffer(depth_buffer) {}

    Matrix<float> vertex(int iface, int nthvert) override{
        Vec3f v = model->vert(iface, nthvert);
        // Set the column of varying_uv to texture position in Vec2f
        varying_uv.set_col(nthvert, model->texcoord(iface, nthvert));
        // Set the column of vertex the triangle using vert index
        Matrix<float> transformed_vert = Viewport*uniform_M*homogonize(v, 1.);
        varying_tri.set_col(nthvert, dehomogonize(transformed_vert));

        return transformed_vert;
    }
    // bar is the barycentric of that vertex
    bool fragment(Vec3f bar, TGAColor &color) override {
        // Convert barycentric vector to a matrix
        // NOTE: Somehow making a new variable is faster than making it inline?
        Matrix<float> bary = Matrix(bar); // 1x3 row matrix that represent a vector
        // Matrix<float> uv = varying_uv*Matrix(bar) <- Slower!
        Matrix<float> mat_uv = varying_uv*bary; // 1x2 Matrix (Basically a Vec2f)

        // Get shadow position from buffer
        Matrix<float> matrix_p = varying_tri * bary;
        Vec3f p = { matrix_p[0][0], matrix_p[1][0], matrix_p[2][0]};
        Matrix<float> shadow_buffer_pt = uniform_Mshadow* homogonize(p, 1.f);
        Vec3f shadow_p = dehomogonize(shadow_buffer_pt);
        auto shadow_buf_idx = static_cast<size_t>(shadow_p.x + shadow_p.y * width);
        float shadow = .3 + 7*(depth_buffer[shadow_buf_idx] < shadow_p.z+43.34); // magic coeff to avoid z-fighting

        Vec2f uv = Vec2f(mat_uv[0][0], mat_uv[1][0]);
        // Get the normal vector of that mesh based on the setting
        Vec3f norm;
        if(!use_normal)
            norm = model->normal(uv.u, uv.v);
        else {
            TGAColor normal_color = normal_file.get(uv.u * normal_file.get_width(), uv.v * normal_file.get_height());
            norm = Vec3f(normal_color.r, normal_color.g, normal_color.b);
        }
        /// Insanely costly calculations
        // Transform the normal vector to the eye space
        Vec3f n = dehomogonize(uniform_MIT*homogonize(norm, 0.f)).normalize();
        // Same as above
        Vec3f l = dehomogonize(uniform_M  *homogonize(light, 0.f)).normalize();
        l.z = -l.z; // I think this formula is meant to work for a different axis?
        l.y = -l.y;
        l.x = -l.x;
        Vec3f r = (n*(n*l*2.f) - l).normalize(); // reflection vector
        float diff = std::max(0.f, n * l); // diffuse intensity value
        // Specular
        float spec = 0.f;
        if(use_specular) {
            float spec_map_val = specular_file.get(uv.u * specular_file.get_width(), uv.v * specular_file.get_height()).r;
            spec = pow(std::max(r.z, 0.f), spec_map_val);
        } else {
            spec = std::max(r.z, 0.f);
        }
        TGAColor texColor = tex_file.get(uv.u * tex_file.get_width(), uv.v * tex_file.get_height());
        for (int i = 0; i < 3; i++) color[i] = std::min<float>(5 + texColor[i]*shadow*(1.2*diff + .6*spec), 255.f);
        color = texColor * diff;

        return false;
    }
};

// The famous Rainbow Triangle
// vertex shader is discarded entirely
struct RainbowShader: IShader {
    Matrix<float> vertex(int iface, int nthvert) override {
        Matrix<float> dummy = Matrix4x4f();
        return dummy;
    }

    bool fragment(Vec3f bar, TGAColor &color) override{
        TGAColor rainbow(bar.x * 255, bar.y * 255, bar.z * 255, 255);
        color = rainbow;
        return false;
    }
};

// Copy zbuffer to a framebuffer (Image in this case)
struct DepthShader: IShader {
    Matrix3x3<float> varying_tri; // 3x3 matrix containing vertex position of a trig

    // Typical vertex rendering
    Matrix<float> vertex(int iface, int nthvert) override{
        Vec3f v = model->vert(iface, nthvert);
        // Set the column of varying_uv to texture position in Vec2f
        Matrix<float> transformed_vert = Viewport*uniform_M*homogonize(v, 1.);
        varying_tri.set_col(nthvert, dehomogonize(transformed_vert));
        return transformed_vert;
    }

    //
    bool fragment(Vec3f bar, TGAColor &color) override {
        Matrix<float> bary(bar);
        Matrix<float> pt = varying_tri * bary; // point interpolated after transformation

        // Set the brightness based on how far is it from the camera
        // clamp
        float dist = std::clamp(pt[2][0]/depth, 0.f, 1.f);
        color = TGAColor(255, 255, 255) *
                (dist);

        return false;
    }
};
#endif //SHADERS_HPP
