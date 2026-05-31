//
// Created by iuv on 7/30/25.
//


#ifndef RENDER_HPP
#define RENDER_HPP
#include <optional>
#include "../../Resources/tgaimage.h"
#include "../../Utilities/Math/geometry.h"
#include "../../Utilities/RectBuffer.h"

class IShader;
class SDL_Surface;
class Model;



class Renderer {
public:
    Renderer() { }
    ~Renderer();
    //void LoadRendererSettings();
    void Render(SDL_Surface *surface);
    void DestroyRenderer();
    void SetCamera(const Vec3f &eye, const Vec3f &up, const Vec3f &cam, float fov);
    void SetupScene(int width, int height, float depth, const Vec3f &light);
    void SetupMaterial(std::string texPath, std::string normalPath, std::string specPath);
    void SetupModel(std::string modelPath);
    void ChangeResolution(int width, int height);
    void SetupBuffers();

    Model *mModel;
    Material *mMaterial;
    Scene *mScene;
    Camera *mCamera;
    RectBuffer *mShadowBuffer;
    RectBuffer *mZbuffer;
private:
    Vec3f rasterize(IShader *shader, int iface, int nthvert);
};
#endif //RENDER_HPP
