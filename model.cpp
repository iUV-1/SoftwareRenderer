//
// Created by iUV on 9/7/2024.
//
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include "model.h"

#include "my_gl.hpp"

Model::Model(const char *filename) {
    std::ifstream in;
    in.open (filename, std::ifstream::in);
    if (in.fail()) return;
    std::string line;
    while (!in.eof()) {
        std::getline(in, line);
        std::istringstream iss(line.c_str());
        char trash;
        if (!line.compare(0, 2, "v ")) {
            iss >> trash;
            Vec3f v;
            for (int i=0;i<3;i++) iss >> v.raw[i];
            verts_.push_back(v);
        } else if (!line.compare(0, 2, "f ")) {
            // v_idx is the vertex index from that face
            // t_idx is the texture index from that face
            // n_idx is the normal index from that face
            std::vector<int> f;
            std::vector<int> t;
            std::vector<int> n;
            int v_idx, t_idx, n_idx;
            iss >> trash;
            while (iss >> v_idx >> trash >> t_idx >> trash >> n_idx) {
                v_idx--; // in wavefront obj all indices start at 1, not zero
                t_idx--;
                n_idx--;
                f.push_back(v_idx);
                t.push_back(t_idx);
                n.push_back(n_idx);
            }
            faces_.push_back(f);
            faces_texture.push_back(t);
            faces_normal.push_back(n);
        } else if (!line.compare(0, 3, "vt ")) {
            // trash vt because trash is a char so we need to do this twice
            iss >> trash;
            iss >> trash;

            Vec2f uv;
            iss >> uv.u >> uv.v;
            texcoords_.push_back(uv);
        } else if (!line.compare(0, 3, "vn ")) {
            iss >> trash;
            iss >> trash;
            Vec3f normal;
            iss >> normal.x >> normal.y >> normal.z;
            normals_.push_back(normal);
        }
    }
    std::cerr << "# v# " << verts_.size() << " f# "  << faces_.size() << std::endl;
}

Model::~Model() {
}

int Model::nverts() {
    return (int)verts_.size();
}

int Model::nfaces() {
    return (int)faces_.size();
}

std::vector<int> Model::face(int idx) {
    return faces_[idx];
}

std::vector<int> Model::face_tex(int idx) {
    return faces_texture[idx];
}

Vec3f Model::vert(int i) {
    return verts_[i];
}

Vec3f Model::vert(int iface, int nthvert) {
    return verts_[faces_[iface][nthvert]];
}

Vec3f Model::normal(int iface, int nthvert) {
    return normals_[faces_[iface][nthvert]];
}

Vec2f Model::texcoord(int i) {
    return texcoords_[i];
}
