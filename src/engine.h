//
// Created by iUV on 5/12/2026.
//

#pragma once
#include "string"
#include "OS/Window.h"

class Camera;
class Scene;
class IRenderer;
//constexpr float depth = 255.f;

class Engine {
private:
    void InitializeUniforms();
    Window *mWindow;
    const bool *keystate;
    const Vec2f *mRelativeMouse;

    Camera *mCamera;
    Scene *mScene;
    IRenderer *mRenderer;

    double mRenderTime;
    double mUpdateTime;
public:
    void Initialize(int agrc, char *argv[]);
    void Update();
    void RenderIMGUI();
    void HandleInput();
};


