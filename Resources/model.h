//
// Created by iUV on 9/7/2024.
//


#pragma once

#include <vector>
#include "../Utilities/Math/geometry.h"

class Model {
private:
    std::vector<Vec3f> verts_;
    std::vector<std::vector<int> > faces_;
    std::vector<std::vector<int>> faces_texture;
    std::vector<std::vector<int>> faces_normal;
    std::vector<Vec2f> texcoords_;
    std::vector<Vec3f> normals_;
public:
    Model(const char *filename);
    ~Model();
    int nverts();
    int nfaces();
    Vec3f vert(int i);
    Vec3f vert(int iface, int nthvert);
    Vec3f normal(int i);
    Vec3f normal(int iface, int nthvert);
    Vec2f texcoord(int i);
    Vec2f texcoord(int iface, int nthvert);
    std::vector<int> face(int idx);
    std::vector<int> face_tex(int idx);
};
