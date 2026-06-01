//
// Created by iUV on 5/12/2026.
//

#include "engine.h"
#include "string"
#include "Utilities/CycleTimer.hpp"
#include "Utilities/Interfaces/IRenderer.h"
#include "Resources/Camera.h"
#include "Resources/Scene.h"
#include <stdexcept>

constexpr int width = 800;
constexpr int height = 800;

using namespace std;

void Engine::InitializeUniforms()
{

}

void Engine::Initialize(int argc, char *argv[])
{
    double setupTimeTaken = CycleTimer::currentSeconds();
    // --- WINDOW SETUP ---
    mWindow = new Window("Software Renderer", width, height);
    keystate = mWindow->GetKeyboardState();
    mRelativeMouse = mWindow->GetRelativeMouse();
    // --- SCENE SETUP ---
    string modelPath = "obj/african_head.obj";
    if(argc >= 2) {
        modelPath = argv[1];
    }
    Mesh *mesh = new Mesh(modelPath.c_str());

    string diffPath = "obj/UV Grid.tga";
    string normalPath;
    string specularPath;

    if(argc >= 3)
        diffPath = argv[2];
    if(argc > 3)
        normalPath = argv[3];
    if(argc > 4)
        specularPath = argv[4];

    Material *material = new Material(diffPath, normalPath, specularPath);

//    Model *model = new Model();
    // --- RENDERER SETUP ---

    // Timing
    setupTimeTaken = CycleTimer::currentSeconds() - setupTimeTaken;
    SDL_Log("Setup time: %.2f ms\n", 1000.f * setupTimeTaken);
}

void Engine::RenderIMGUI()
{
    ImGui::Begin("Renderer settings");

    ImGui::Text("Render time/FPS: %.3f ms/frame (%.1f FPS)", mRenderTime * 1000, 1/mRenderTime);
    ImGui::Text("Update time/FPS: %.3f ms/frame (%.1f FPS)", mUpdateTime * 1000, 1/mUpdateTime);

    // Since Light is a continuous array, &rScene->Light.x works
    //ImGui::DragFloat3("Light", &rScene->Light.x, 0.1, -1.0, 1.0, "%.3f");
    ImGui::End();
}

void Engine::HandleInput()
{
    Vec2i result;
    result.x = keystate[SDL_SCANCODE_A] - keystate[SDL_SCANCODE_D];
    result.y = keystate[SDL_SCANCODE_W] - keystate[SDL_SCANCODE_S];

    Vec2f mouse = *mRelativeMouse;
    mCamera->Eye.x += 0.1f * result.x;
    mCamera->Eye.z += 0.1f * result.y;
    mCamera->Cam.x += 0.1f * result.x - mouse.x * 0.001;
    mCamera->Cam.z += 0.1f * result.y;
    mCamera->Cam.y -= mouse.y * 0.001;
}

void Engine::Update()
{
    mUpdateTime = CycleTimer::currentSeconds();

    HandleInput();
    // --- RENDER ---
    mRenderTime = CycleTimer::currentSeconds();
    vector<Model*> models = mScene->GetModels();
    for(int i = 0; i < mScene->GetModelSize(); i++) {
        mRenderer->RenderModel(models[i]);
    }
    mRenderTime = CycleTimer::currentSeconds() - mRenderTime;


    RenderIMGUI();
    mWindow->Update();
}
