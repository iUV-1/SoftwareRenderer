//
// Created by iUV on 5/12/2026.
//

#include "engine.h"
#include "string"
#include "Utilities/CycleTimer.hpp"
#include "Utilities/Interfaces/IRenderer.h"
#include "Resources/Camera.h"
#include "Resources/Scene.h"
#include "Utilities/Math/math.h"
#include "Renderer/SoftwareRenderer/render.hpp"
#include <stdexcept>

constexpr int window_width = 800;
constexpr int window_height = 800;

using namespace std;

void Engine::InitializeUniforms()
{
    Matrix4x4f projection = Project(60, mWindow->mWidth, mWindow->mHeight, 0.1f, 255.f);
    Matrix4x4f viewport = SetViewport(mWindow->mWidth / 8, mWindow->mHeight / 8, mWindow->mWidth * 3./4, mWindow->mHeight * 3./4, 255.f);
    Matrix4x4f modelview = LookAt(mCamera->Eye, mCamera->Cam, mCamera->Up);
    dynamic_cast<Renderer*>(mRenderer)->SetupUniforms(projection, modelview, viewport);
}

void Engine::ResizeWindow(int width, int height)
{
    mCamera->Width = width;
    mCamera->Height = height;
    dynamic_cast<Renderer*>(mRenderer)->SetupBuffers(mWindow->mWidth, mWindow->mHeight);
}

void Engine::Initialize(int argc, char *argv[])
{
    double setupTimeTaken = CycleTimer::currentSeconds();
    // --- WINDOW SETUP ---
    mWindow = new Window("Software Renderer", window_width, window_height);
    keystate = mWindow->GetKeyboardState();
    mRelativeMouse = mWindow->GetRelativeMouse();
    // --- SCENE SETUP ---
    string modelPath = "../obj/african_head.obj";
    if(argc >= 2) {
        modelPath = argv[1];
    }
    Mesh mesh(modelPath.c_str());

    string diffPath = "../obj/UV Grid.tga";
    string normalPath;
    string specularPath;

    if(argc >= 3)
        diffPath = argv[2];
    if(argc >= 4)
        normalPath = argv[3];
    if(argc >= 5)
        specularPath = argv[4];

    Material material(diffPath, normalPath, specularPath);
    mCamera = new Camera(
            Vec3f(1, 1, 3),
            Vec3f(0, 1, 0),
            Vec3f(-1, 1, -1),
            60
            );

    Model* model = new Model(Vec3f(0, 0, 0), mesh, material);
    mScene = new Scene();
    mScene->AddModel( model );
    mScene->SetCamera(mCamera);
    mScene->SetLight(Vec3f(1, -1, 0));
    // --- RENDERER SETUP ---
    mRenderer = new Renderer(mWindow);
    dynamic_cast<Renderer*>(mRenderer)->SetupBuffers(mWindow->mWidth, mWindow->mHeight);
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
    Vec3f *light = mScene->GetLight();
    ImGui::DragFloat3("Light", &light->x, 0.1, -1.0, 1.0, "%.3f");
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
    // Capture the start time in a temporary local variable
    // This is done for mUpdateTime since we render imgui in the middle of this block
    // meaning when imgui renders, it will be currentSeconds(), which is the amount of seconds have passed since this program starts.
    // imgui wont get to render the proper updateTime since the calculation is after RenderImGUI()
    double startTime = CycleTimer::currentSeconds();

    mWindow->StartUpdate();
    HandleInput();
    InitializeUniforms();
    // --- RENDER ---
    mRenderTime = CycleTimer::currentSeconds();
    //dynamic_cast<Renderer*>(mRenderer)->SetupBuffers(mWindow->mWidth, mWindow->mHeight);
    mRenderer->RenderScene(mScene);
    mRenderTime = CycleTimer::currentSeconds() - mRenderTime;

    RenderIMGUI();
    mWindow->EndUpdate();
    mUpdateTime = CycleTimer::currentSeconds() - startTime;
}

Engine::~Engine() {
    delete mWindow;
    delete mRenderer;
    delete mScene;
}