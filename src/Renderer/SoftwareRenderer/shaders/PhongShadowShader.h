//
// Created by iUV on 1/8/2026.
//
#pragma once
#include <iostream>
#include <random>
#include <algorithm>

#include "../../../Resources/Mesh.h"
#include "../../../Utilities/Interfaces/IShader.hpp"
#include "../../../Utilities/Math/geometry.h"
#include "../my_gl.hpp"
#include "../render.hpp"

// Like above but includes a shadow pass
struct PhongShaderShadow: IShader {
    Matrix3x3<float> varying_tri; // 3x3 matrix containing verticies of a trig
    Matrix3x3<float> varying_shadow_tri; // 3x3 matrix containing verticies of a shadow trig
    Matrix4x4f uniform_Mshadow; // Shadow transformation
    Vec3f uniform_light;
    RectBuffer &depth_buffer; // for shadow
    Matrix4x4f Viewport;
    Matrix4x4f uniform_M;
    Matrix<float, 3, 2> varying_uv;
    Matrix4x4f uniform_MIT; // Invert transpose

    Mesh *mesh;
    Material *mat;
    PhongShaderShadow(Matrix4x4f viewport, Matrix4x4f projection, Matrix4x4f modelview,
                      Model *model, Vec3f light,
                      Matrix4x4f uniform_shadow, RectBuffer &depth_buffer)
    : Viewport(viewport), depth_buffer(depth_buffer), uniform_light(light) {
        uniform_M = projection*modelview;
        uniform_MIT = uniform_M;
        uniform_MIT.inverseTranspose();
        Matrix4x4f M = (viewport*projection*modelview);
        M.invert();
        uniform_Mshadow = uniform_shadow*M;
        mesh = &model->mMesh;
        mat = &model->mMaterial;
    }

    Vec4f vertex(int iface, int nthvert) override{
        Vec3f v = mesh->vert(iface, nthvert);
        // Set the column of varying_uv to texture position in Vec2f
        varying_uv.set_col(nthvert, mesh->texcoord(iface, nthvert));
        // Set the column of vertex the triangle using vert index
        Vec4f transformed_vert = Viewport*uniform_M*Vec4f(v, 1.);
//        bool x_c = -transformed_vert.w <= transformed_vert.x && transformed_vert.x <= transformed_vert.w;
//        bool y_c = -transformed_vert.w <= transformed_vert.y && transformed_vert.y <= transformed_vert.w;
//        bool z_c = 0 <= transformed_vert.z && transformed_vert.z <= transformed_vert.w;
//        if(!x_c || !y_c || !z_c) return Vec4f(1, 1, 1, 1);
        varying_tri.set_col(nthvert, Vec3f(transformed_vert));
        varying_shadow_tri.set_col(nthvert, Vec3f(uniform_Mshadow * Vec4f(v, 1.)));

        return transformed_vert;
    }
    // bar is the barycentric of that vertex
    bool fragment(Vec3f bar, TGAColor &color) override {
        // Interpolate the uv vector
        // Heavy operation, Matrix*Vec<3> returns a matrix
        // And then this cast it to a vector
//        Vec2f uv(varying_uv * bar);

        Vec2f uv{};
        uv.x = varying_uv(0,0) * bar.x + varying_uv(1, 0) * bar.y + varying_uv(2, 0) * bar.z;
        uv.y = varying_uv(0, 1) * bar.x + varying_uv(1, 1) * bar.y + varying_uv(2, 1) * bar.z;
        // Get shadow position from buffer
        Vec3f p = varying_tri * bar;
        Vec3f shadow_p = Vec3f(uniform_Mshadow * Vec4f(p, 1.));
        auto shadow_buf_idx = int(shadow_p.x)  + int(shadow_p.y) * depth_buffer.width;

        // Get the normal vector of that mesh based on the setting
        Vec3f norm;
        if(!mat->UseNormal)
            norm = mesh->normal(uv.u, uv.v);
        else {
            TGAColor normal_color = TEX2D(mat->NormalFile, uv);
            norm = Vec3f(normal_color.r, normal_color.g, normal_color.b);
        }
        /// Insanely costly calculations
        // Transform the normal vector to the mEye space
        Vec3f n = Vec3f(uniform_MIT*Vec4f(norm, 1.f)).normalize();
        // Same as above
        Vec3f l = Vec3f(uniform_M  *Vec4f(uniform_light, 1.f)).normalize();
        l *= -1;// The reflection formula below is for object pointing to the light.
        // Shadow
        float slope_bias = std::max( 0.5f* (1.0f - n * (l*-1)), 1.f );
        //slope_bias = 43.34f;
        float depth_p = depth_buffer[shadow_buf_idx];
        if(depth_p != -MAX_FLOAT) {
            depth_p -= slope_bias;
        }

        float shadow = (depth_p < shadow_p.z) ? 1.f : 0.3f;

        Vec3f r = (n*(n*l*2.f) - l).normalize(); // reflection vector
        float diff = std::max(0.f, n * l); // diffuse intensity value
        // Specular
        float spec = 0.f;
        if(mat->UseSpecular) {
            float spec_map_val = TEX2D(mat->SpecularFile, uv).r;
            spec = std::pow(std::max(r.z, 0.f), spec_map_val);
        } else {
            spec = std::max(r.z, 0.f);
        }

        TGAColor texColor = TEX2D(mat->TexFile, uv);
        for (int i = 0; i < 3; i++)
            color[i] = std::min<float>(5 + texColor[i]*shadow*(1.2*diff + .6*spec), 255.f);

        return false;
    }
};
