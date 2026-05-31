//
// Created by iUV on 5/27/2026.
//

#ifndef SOFTWARERENDERER_CAMERA_H
#define SOFTWARERENDERER_CAMERA_H
#include "../Utilities/Math/geometry.h"

class Camera {
public:
    Camera() { }
    // Copy constructor
    //Camera(const Camera *camera);
    Camera(const Vec3f &eye, const Vec3f &up, const Vec3f &cam, float fov): Eye(eye), Up(up), Cam(cam), Fov(fov) {}
    //~Camera();
    Vec3f Eye;
    Vec3f Up;
    Vec3f Cam;
    float Depth;

    float Fov;
};

#endif //SOFTWARERENDERER_CAMERA_H
