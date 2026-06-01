//
// Created by iUV on 5/27/2026.
//

#ifndef SOFTWARERENDERER_SCENE_H
#define SOFTWARERENDERER_SCENE_H
#include "../../Utilities/Math/geometry.h"
#include "Model.h"
#include "Camera.h"

class Scene {
private:
    std::vector<Model*> mModels;
    Vec3f mLight;
    Camera *mCamera;
public:
    Scene() { }
    Scene(const Scene *scene);
//    Scene(int width, int height, float depth, const Vec3f &light)
//            : Width(width), Height(height), Light(light) {}
    //~Scene();
//    std::vector<Vec3f> Light; // Only supports directional light for now
    void AddModel(Model *model);
    void SetCamera(Camera *camera);
    std::vector<Model*> GetModels() {
        return mModels;
    }
    int GetModelSize() {
        return mModels.size();
    };
};


#endif //SOFTWARERENDERER_SCENE_H
