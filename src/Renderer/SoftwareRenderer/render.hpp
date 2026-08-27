//
// Created by iuv on 7/30/25.
//


#ifndef RENDER_HPP
#define RENDER_HPP
#include <optional>
#include "../../Resources/tgaimage.h"
#include "../../../Utilities/Math/geometry.h"
#include "../../../Utilities/RectBuffer.h"
#include "../../../Utilities/Interfaces/IRenderer.h"

class IShader;
class SDL_Surface;
class Mesh;
class Scene;
class Camera;

class Renderer: public IRenderer {
public:
    //Renderer() { }
    Renderer(Window *window):  mWindow(window) {
        SetupBuffers(mWindow->mWidth, mWindow->mHeight);
    };
    ~Renderer();
    //void LoadRendererSettings();
    void Render(SDL_Surface *surface);
    void SetupBuffers(int width, int height);
    void RenderScene(Scene *scene) override;
    void RenderIMGUI() override;
    void SetupUniforms(Matrix4x4f const &projection, Matrix4x4f const &modelview, Matrix4x4f const &viewport);
    Scene *mScene;
    Camera *mCamera;
    Window *mWindow;
    RectBuffer *mShadowBuffer = nullptr;
    RectBuffer *mZbuffer = nullptr;

    // Own by the engine
    Matrix4x4f mProjection;
    Matrix4x4f mModelView;
    Matrix4x4f mViewport;
private:
    Vec3f rasterize(IShader *shader, int iface, int nthvert);
};
#endif //RENDER_HPP
