//
// Created by iUV on 5/26/2026.
//

#include "math.h"
Matrix4x4f ModelView;
Matrix4x4f Projection;
Matrix4x4f Viewport;

inline float degToRad(float deg) {
    return deg * M_PI / 180;
}

// Right hand side
void Project(float fov, float width, float height, float near, float far) {

    double fovCoeff = 1. / tan(degToRad(fov) / 2); // cot(fov/2)
    float nearFar = near - far;
    float aspectRatio = width/height;
    Projection[0][0] = 1/aspectRatio * fovCoeff;
    Projection[1][1] = fovCoeff;
    Projection[2][2] = far/nearFar;
    Projection[2][3] = -1;
    Projection[3][2] = -far * near/nearFar;
    Projection[3][3] = 0;
}

/// Set the viewport plane matrix
/// x: start point in x
/// y: start point in y
/// w: width
/// h: height
/// depth: depth
void SetViewport(int x, int y, float w, float h, float depth) {
    Viewport = Matrix4x4f::identity();
    Viewport[0][3] = x + w / 2.f;
    Viewport[1][3] = y + h / 2.f;
    Viewport[2][3] = depth / 2.f;

    Viewport[0][0] = w / 2.f;
    Viewport[1][1] = h / 2.f;
    Viewport[2][2] = depth / 2.f;
}

// Similar to gluLookAt, create a camera transformation matrix
// Formula (8.4) in textbook
void LookAt(Vec3f eye, Vec3f center, Vec3f up) {
    Vec3f z = (eye-center).normalize();
    Vec3f x = (up^z).normalize();
    Vec3f y = (z^x).normalize();
    // Inverse rotation matrix
    auto Minv = Matrix4x4f::identity();
    // Translation matrix
    auto Tr = Matrix4x4f::identity();
    for (int i = 0; i < 3; i++) {
        Minv[0][i] = x[i];
        Minv[1][i] = y[i];
        Minv[2][i] = z[i];
        Tr[i][3] = -eye[i];
    }
    ModelView = Minv * Tr;
}

void OrthoProject(float coeff) {
    Projection = Matrix4x4f ::identity();
    Projection[3][2] = 0;
}