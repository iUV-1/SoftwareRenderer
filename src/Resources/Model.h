//
// Created by iUV on 5/31/2026.
//

#ifndef SOFTWARERENDERER_MODEL_H
#define SOFTWARERENDERER_MODEL_H

#include "Mesh.h"
#include "Material.h"

class Model {
private:


public:
    Vec3f mPosition;
    Mesh mMesh;
    Material mMaterial;
    Model() = delete;
    Model(const Vec3f &pos, const Mesh &mesh, const Material &mat):
         mPosition(pos), mMesh(mesh), mMaterial(mat) { }
};


#endif //SOFTWARERENDERER_MODEL_H
