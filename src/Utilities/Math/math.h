//
// Created by iUV on 5/26/2026.
//

#pragma once

#include "geometry.h"

Matrix4x4f SetViewport(int x, int y, float w, float h, float depth);
Matrix4x4f LookAt(Vec3f eye, Vec3f center, Vec3f up);
Matrix4x4f OrthoProject(float coeff);
Matrix4x4f Project(float fov, float width, float height, float near, float far);