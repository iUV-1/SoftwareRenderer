//
// Created by iUV on 5/26/2026.
//
// For MSVC because they don't define M_PI by default
#define _USE_MATH_DEFINES
#include "math.h"
//Matrix4x4f Viewport;

inline float degToRad(float deg) {
    return deg * M_PI / 180;
}

// Right hand side
Matrix4x4f Project(float fov, float width, float height, float near, float far) {
    Matrix4x4f Projection;
    double fovCoeff = 1. / tan(degToRad(fov) / 2); // cot(fov/2)
    float nearFar = near - far;
    float aspectRatio = width/height;
    Projection(0,0) = 1/aspectRatio * fovCoeff;
    Projection(1,1) = fovCoeff;
    Projection(2,2) = far/nearFar;
    Projection(3, 2) = -1;
    Projection(2, 3) = -far * near/nearFar;
    Projection(3, 3) = 0;
    return Projection;
}

/// Set the viewport plane matrix
/// x: start point in x
/// y: start point in y
/// w: width
/// h: height
/// depth: depth
Matrix4x4f SetViewport(int x, int y, float w, float h, float depth) {
    Matrix4x4f Viewport;
    Viewport = Matrix4x4f::identity();
    Viewport(3, 0) = x + w / 2.f;
    Viewport(3, 1) = y + h / 2.f;
    Viewport(3, 2) = depth / 2.f;

    Viewport(0, 0) = w / 2.f;
    Viewport(1, 1) = h / 2.f;
    Viewport(2, 2) = depth / 2.f;
    return Viewport;
}

// Similar to gluLookAt, create a camera transformation matrix
// Formula (8.4) in textbook
Matrix4x4f LookAt(Vec3f eye, Vec3f center, Vec3f up) {
    Matrix4x4f ModelView;
    Vec3f z = (eye-center).normalize();
    Vec3f x = (up^z).normalize();
    Vec3f y = (z^x).normalize();
    // Inverse rotation matrix
    auto Minv = Matrix4x4f::identity();
    // Translation matrix
    auto Tr = Matrix4x4f::identity();
    for (int i = 0; i < 3; i++) {
        Minv(i, 0) = x[i];
        Minv(i ,1) = y[i];
        Minv(i, 2) = z[i];
        Tr(3, i) = -eye[i];
    }
    ModelView = Minv * Tr;
    return ModelView;
}

Matrix4x4f OrthoProject(float coeff) {
    Matrix4x4f Projection;
    Projection = Matrix4x4f ::identity();
    Projection(2, 3) = coeff;
    return Projection;
}