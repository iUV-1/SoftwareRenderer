//
// Created by iUV on 1/25/2026.
//

#pragma once
#include <iostream>
#include <random>
#include <algorithm>

#include "../model.h"
#include "../shaders.hpp"
#include "../geometry.h"
#include "../my_gl.hpp"
#include "../render.hpp"
// Gouraud Shader uses vertex data to calculate light value
// However, this shader uses interpolation to compute the normal vector per pixel
// This is also called smooth shading
//struct GouraudShader: IShader {
//    Vec3f varying_intensity; // intensity of a vertex
//    Matrix<float> varying_uv = Matrix<float>(2, 3); // 2x3 matrix containing uv coordinate of 3 vertex (a trig)
//    Matrix4x4f uniform_M; // Projection*ModelView
//    Matrix4x4f uniform_MIT; // same as above but invert_transpose()
//
//    Vec4f vertex(int iface, int nthvert) override{
//        //light.normalize();
//        Vec3f v = model->vert(iface, nthvert);
//        Vec3f n = model->normal(iface, nthvert);
//        // Set the column of varying_uv to texture position in Vec2f
//        varying_uv[0][nthvert] = model->texcoord(iface, nthvert).x;
//        varying_uv[1][nthvert] = model->texcoord(iface, nthvert).y;
//        // Cap at 0
//        varying_intensity[nthvert] = std::max(0.f, n*light);
//        return Viewport*Projection*ModelView*homogonize(v, 1.);
//    }
//    // bar is the barycentric of that vertex
//    bool fragment(Vec3f bar, TGAColor &color) override {
//        // Cap at 0
//        // Convert barycentric vector to a matrix
//        // NOTE: Somehow making a new variable is faster than making it inline?
//        Matrix<float> bary = Matrix(bar); // 1x3 row matrix that represent a vector
//        // Matrix<float> uv = varying_uv*Matrix(bar) <- Slower!
//        Matrix<float> uv = varying_uv*bary; // 1x2 Matrix (Basically a Vec2f)
//
//        Vec3f norm;
//        if(false)
//            norm = model->normal(uv[0][0], uv[1][0]);
//        else {
//            TGAColor normal_color = normal_file.get(uv[0][0] * normal_file.get_width(), uv[1][0] * normal_file.get_height());
//            norm = Vec3f(normal_color.r, normal_color.g, normal_color.b);
//        }
//        // Transforming
//        Vec3f l = dehomogonize(uniform_M  *homogonize(light, 0.f)).normalize();
//        //l.z = -l.z;
//        l.x = -l.x;
//        l.y = -l.y;
//        l.z = -l.z;
//        Vec3f n = dehomogonize(uniform_MIT*homogonize(norm, 0.f)).normalize();
//        if(varying_intensity*bar != n*l) {
//            std::cout << "mismatched value" << std::endl;
//        }
//
//        float diff = std::max(0.f, n * l); // diffuse intensity value
//
//        TGAColor texColor = mTexFile.get(uv[0][0] * mTexFile.get_width(), uv[1][0] * mTexFile.get_height());
//        color = texColor * diff;
//        return false;
//    }
//};