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

using namespace std;

void Engine::InitializeUniforms()
{

}

void Engine::Initialize(int agrc, char *argv[])
{
    double beforeSetupTime = CycleTimer::currentSeconds();
    // --- WINDOW SETUP ---
    mWindow = new Window("Software Renderer", 800, 800);
    keystate = mWindow->GetKeyboardState();
    mRelativeMouse = mWindow->GetRelativeMouse();
    // --- SCENE SETUP ---

    // --- RENDERER SETUP ---

    // Timing
    double setupTimeTaken = CycleTimer::currentSeconds() - beforeSetupTime;
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
