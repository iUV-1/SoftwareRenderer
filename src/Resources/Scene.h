//
// Created by iUV on 5/27/2026.
//

#ifndef SOFTWARERENDERER_SCENE_H
#define SOFTWARERENDERER_SCENE_H
#include <memory>
#include "../../Utilities/Math/geometry.h"
#include "Model.h"
#include "Camera.h"


class Scene {
private:
    std::vector< Model* > mModels;
    Vec3f mLight;
    Camera *mCamera;
public:
    Scene() { }
    Scene(const Scene *scene);
//    Scene(int width, int height, float depth, const Vec3f &light)
//            : Width(width), Height(height), Light(light) {}
    //~Scene();
//    std::vector<Vec3f> Light; // Only supports directional light for now

// Now scene holds this model
    void AddModel(Model *model) {
        mModels.push_back( model );
    };

    void SetCamera(Camera *camera) {
        mCamera = camera;
    };

    void SetLight(Vec3f light) {
        mLight = light;
    }
    std::vector< Model* >& GetModels() {
        return mModels;
    }
    int GetModelSize() {
        return mModels.size();
    };

    Vec3f* GetLight() {
        return &mLight;
    }

    Camera* GetCamera() {
        return mCamera;
    }
};


#endif //SOFTWARERENDERER_SCENE_H
