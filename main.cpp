//
// Created by iUV on 1/6/2026.
//
#define SDL_MAIN_USE_CALLBACKS 1

#include "pch.h"
#include <SDL3/SDL_main.h>
#include <imgui.h>
#include <imgui_impl_sdlrenderer3.h>
#include <imgui_impl_sdl3.h>
#include <string>

#include "render.hpp"
#include "CycleTimer.hpp"

SDL_Window *window;
//SDL_Surface *windowSurface;
SDL_Renderer *sdlRenderer;

int width = 800;
int height = 600;
constexpr float depth = 255.f;
Renderer *renderer;

struct GlobalAppState {
    const bool *keystate;
    ImGuiIO imGuiIo;
};

SDL_AppResult SDL_AppInit(void **appstate, int argc, char *argv[]) {
    // Get the global struct by casting
    appstate[0] = new GlobalAppState;
    if(!SDL_Init(SDL_INIT_VIDEO | SDL_INIT_JOYSTICK)) {
        return SDL_APP_FAILURE;
    }
    window = SDL_CreateWindow("Software Renderer", width, height, SDL_WINDOW_RESIZABLE);
    // --- SETUP RENDERER ---
    double beforeSetup = CycleTimer::currentSeconds();
    renderer = new Renderer();
    renderer->SetCamera(Vec3f(1, 1, 3),
                         Vec3f(0, 1, 0),
                        Vec3f(0, 0, 0),
                        90);
    if(argc < 2)
        renderer->SetupModel("obj/african_head.obj");
    else
        renderer->SetupModel(argv[1]);

    std::string diffPath;
    std::string normalPath;
    std::string specularPath;
    if(argc < 3)
        diffPath = "obj/UV Grid.tga";
    else
        diffPath = argv[2];

    if(argc >= 3)
        normalPath = argv[3];
    if(argc >= 4)
        specularPath = argv[4];

    renderer->SetupMaterial(diffPath, normalPath, specularPath);
    renderer->SetupBuffers();
    renderer->SetupScene(width, height, depth, Vec3f(1.0, 1.0, 1.0));
    double setupTimeTaken = CycleTimer::currentSeconds() - beforeSetup;
    SDL_Log("Setup time: %.2f ms\n", 1000.f * setupTimeTaken);

    // --- SETUP IMGUI ---
    sdlRenderer = SDL_CreateRenderer(window, nullptr);
    if (renderer == nullptr)
    {
        SDL_Log("Error: SDL_CreateRenderer(): %s\n", SDL_GetError());
        return SDL_APP_FAILURE;
    }
    // imgui stuff
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();    (void)io;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
    ((GlobalAppState*)appstate[0])->imGuiIo = io; // may lord forgives me
    // Setup Dear ImGui style
    ImGui::StyleColorsDark();
    //ImGui::StyleColorsLight();
    float main_scale = SDL_GetDisplayContentScale(SDL_GetPrimaryDisplay());

    // Setup scaling
    ImGuiStyle& style = ImGui::GetStyle();
    style.ScaleAllSizes(main_scale);        // Bake a fixed style scale. (until we have a solution for dynamic style scaling, changing this requires resetting Style + calling this again)
    style.FontScaleDpi = main_scale;        // Set initial font scale. (using io.ConfigDpiScaleFonts=true makes this unnecessary. We leave both here for documentation purpose)
    // Setup Platform/Renderer backends
    ImGui_ImplSDL3_InitForSDLRenderer(window, sdlRenderer);
    ImGui_ImplSDLRenderer3_Init(sdlRenderer);

    const bool* keystate = SDL_GetKeyboardState(nullptr);
    ((GlobalAppState*)appstate[0])->keystate = keystate;

    return SDL_APP_CONTINUE;
}

double renderTime = 30.f;
double iterateTime = 30.f;
void HandleInput(const bool *keystate, Camera *cam) {
    Vec2i result;
    result.x = keystate[SDL_SCANCODE_A] - keystate[SDL_SCANCODE_D];
    result.y = keystate[SDL_SCANCODE_W] - keystate[SDL_SCANCODE_S];

    cam->Eye.x += 0.1f * result.x;
    cam->Eye.z += 0.1f * result.y;
    cam->Cam.x += 0.1f * result.x;
    cam->Cam.z += 0.1f * result.y;
}

SDL_AppResult SDL_AppIterate(void *appstate) {
    double beforeIterate = CycleTimer::currentSeconds();
    // --- IMGUI STUFF ---
    ImGuiIO& io = ((GlobalAppState*)appstate)->imGuiIo; // this work
    const bool* keystate = ((GlobalAppState*)appstate)->keystate;
    // new frame blah blah blah
    ImGui_ImplSDLRenderer3_NewFrame();
    ImGui_ImplSDL3_NewFrame();
    ImGui::NewFrame();
    Scene *rScene = renderer->mScene;
    // 2. Show a simple window that we create ourselves. We use a Begin/End pair to create a named window.
    {
        ImGui::Begin("Renderer settings");

        ImGui::Text("Render time/FPS: %.3f ms/frame (%.1f FPS)", renderTime * 1000, 1/renderTime);
        ImGui::Text("Iterate time/FPS: %.3f ms/frame (%.1f FPS)", iterateTime * 1000, 1/iterateTime);

        // Since Light is a continuous array, &rScene->Light.x works
        ImGui::DragFloat3("Light", &rScene->Light.x, 0.1, -1.0, 1.0, "%.3f");
        ImGui::End();
    }
    renderer->ChangeResolution(width, height);
    //rScene->Light.normalize();
    // IMGui render
    ImGui::Render();
    SDL_SetRenderScale(sdlRenderer, io.DisplayFramebufferScale.x, io.DisplayFramebufferScale.y);
    SDL_RenderClear(sdlRenderer);

    // --- MY RENDERER ---
    double before = CycleTimer::currentSeconds();
    // Get the window surface
    SDL_Surface *windowSurface = SDL_GetWindowSurface(window);
    // Clear the surface
    SDL_ClearSurface(windowSurface, 0xFF, 0xFF, 0xFF, 0xFF);
    // Handle movement input
    Camera *rCamera = renderer->mCamera;
    HandleInput(keystate, rCamera);
    // Render
    renderer->Render(windowSurface);
    // Flip the surface
    //SDL_FlipSurface(windowSurface, SDL_FLIP_VERTICAL);
    // Get the timing and FPS
    renderTime = CycleTimer::currentSeconds() - before;
    // Create a texture of the finished result
    SDL_Texture *cpuSurface = SDL_CreateTextureFromSurface(sdlRenderer, windowSurface);
    // Render the texture
    SDL_RenderTexture(sdlRenderer, cpuSurface, nullptr, nullptr);
    ImGui_ImplSDLRenderer3_RenderDrawData(ImGui::GetDrawData(), sdlRenderer);

    SDL_RenderPresent(sdlRenderer);
    iterateTime = CycleTimer::currentSeconds() - beforeIterate;
    return SDL_APP_CONTINUE;
}

SDL_AppResult SDL_AppEvent(void *appstate, SDL_Event *event) {
    ImGui_ImplSDL3_ProcessEvent(event); // imgui event handler

    switch(event->type) {
        case SDL_EVENT_QUIT:
            return SDL_APP_FAILURE; // doing this will quit the program
            break;
        case SDL_EVENT_KEY_DOWN:
            SDL_Log("%d was pressed", event->key.scancode);
            break;
        case SDL_EVENT_WINDOW_RESIZED:
            width = event->window.data1;
            height = event->window.data2;
            break;
        default:
            break;
    }
    return SDL_APP_CONTINUE;
}

void SDL_AppQuit(void *appstate, SDL_AppResult result) {
    delete (GlobalAppState*)appstate;
    ImGui_ImplSDLRenderer3_Shutdown();
    ImGui_ImplSDL3_Shutdown();
    ImGui::DestroyContext();
    delete renderer;
    SDL_DestroyRenderer(sdlRenderer);
    SDL_DestroyWindow(window);
    //SDL_DestroySurface(windowSurface);
    SDL_Quit();
}
