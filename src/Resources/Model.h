//
// Created by iUV on 5/31/2026.
//

#ifndef SOFTWARERENDERER_MODEL_H
#define SOFTWARERENDERER_MODEL_H

#include "Mesh.h"
#include "Material.h"

class Model {
private:
    Vec3f mPosition;
    Mesh mMesh;
    Material mMaterial;
};


#endif //SOFTWARERENDERER_MODEL_H
