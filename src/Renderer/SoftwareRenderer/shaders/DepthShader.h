//
// Created by iUV on 1/10/2026.
//

#pragma once
//#include <iostream>
//#include <random>
//#include <algorithm>
//
//#include "../model.h"
//#include "../shaders.hpp"
//#include "../geometry.h"
//#include "../my_gl.hpp"
//#include "../render.hpp"
//
//struct DepthShader: IShader {
//    // Typical vertex rendering
//    Vec4f vertex(int iface, int nthvert) override{
//        Vec3f v = input->model.vert(iface, nthvert);
//        // Set the column of varying_uv to texture position in Vec2f
//        // Matrix<float> transformed_vert = Viewport*uniform_M*homogonize(v, 1.);
//        return Viewport*uniform_M*homogonize(v, 1.);
//    }
//
//    // Discard the fragment because we are only interested in the zbuffer produced by triangle()
//    bool fragment(Vec3f bar, TGAColor &color) override {
//        color = TGAColor(0,0,0,0);
//        return false;
//    }
//};