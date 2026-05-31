//
// Created by iUV on 1/10/2026.
//

#pragma once
#include <iostream>
#include <random>
#include <algorithm>

#include "../../../Resources/model.h"
#include "../../../Utilities/Interfaces/IShader.hpp"
#include "../../../Utilities/Math/geometry.h"
#include "../my_gl.hpp"
#include "../SoftwareRenderer/render.hpp"

// Copy zbuffer to a framebuffer (Image in this case)
struct DepthShaderImage: IShader {
    Matrix3x3<float> varying_tri; // 3x3 matrix containing vertex position of a trig
    // Typical vertex rendering
    DepthShaderImage(Model *model, Scene *scene): IShader(scene) { uniform_Model = model; uniform_Scene = scene; }
    Vec4f vertex(int iface, int nthvert) override{
        Vec3f v = uniform_Model->vert(iface, nthvert);
        // Set the column of varying_uv to texture position in Vec2f
        Vec4f transformed_vert = Viewport*uniform_M*homogonize(v, 1.);
        varying_tri.set_col(nthvert, dehomogonize(transformed_vert));
        return transformed_vert;
    }

    bool fragment(Vec3f bar, TGAColor &color) override {
        // Set the brightness based on how far is it from the camera
        Vec3f p = varying_tri*bar;
        // clamp
        float dist = std::clamp(p.z/uniform_Scene->Depth, 0.f, 1.f);
        color = TGAColor(255, 255, 255) *
                (dist);

        return false;
    }
};
