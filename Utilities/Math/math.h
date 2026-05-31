//
// Created by iUV on 5/26/2026.
//

#ifndef SOFTWARERENDERER_MATH_H
#define SOFTWARERENDERER_MATH_H

#include "geometry.h"

void SetViewport(int x, int y, float w, float h, float depth);
void LookAt(Vec3f eye, Vec3f center, Vec3f up);
void OrthoProject(float coeff);
void Project(float fov, float width, float height, float near, float far);


#endif //SOFTWARERENDERER_MATH_H
