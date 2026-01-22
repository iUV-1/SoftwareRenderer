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
SDL_Surface *windowSurface;
SDL_Renderer *sdlRenderer;

constexpr int width = 800;
constexpr int height = 800;
constexpr float depth = 255.f;
Renderer *renderer;

Vec2i input;
Vec3f eye(1, 1, 3);
Vec3f up(0, 1, 0);
Vec3f cam(0, 0, 0);

SDL_AppResult SDL_AppInit(void **appstate, int argc, char *argv[]) {
    // Get the global struct by casting
    appstate[0] = new ImGuiIO;
    if(!SDL_Init(SDL_INIT_VIDEO | SDL_INIT_JOYSTICK)) {
        return SDL_APP_FAILURE;
    }
    window = SDL_CreateWindow("Software Renderer", width, height, SDL_WINDOW_RESIZABLE);
    // --- SETUP RENDERER ---
    double beforeSetup = CycleTimer::currentSeconds();
    renderer = new Renderer();

    renderer->SetCamera(eye, up, cam);
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

    Vec3f light(1.0, 1.0, 1.0);
    light.normalize();
    renderer->SetupScene(width, height, 255.f, light);
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
    *(ImGuiIO*)appstate[0] = io; // ig this works
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

    return SDL_APP_CONTINUE;
}

double renderTime = 30.f;
double iterateTime = 30.f;
SDL_AppResult SDL_AppIterate(void *appstate) {
    double beforeIterate = CycleTimer::currentSeconds();
    // --- IMGUI STUFF ---
    ImGuiIO& io = *(ImGuiIO*)appstate; // this work
    // new frame blah blah blah
    ImGui_ImplSDLRenderer3_NewFrame();
    ImGui_ImplSDL3_NewFrame();
    ImGui::NewFrame();
    // 2. Show a simple window that we create ourselves. We use a Begin/End pair to create a named window.
    {
        static float f = 0.0f;
        static int counter = 0;

        ImGui::Begin("Stats");                          // Create a window called "Hello, world!" and append into it.

        ImGui::Text("Render time/FPS: %.3f ms/frame (%.1f FPS)", renderTime * 1000, 1/renderTime);
        ImGui::Text("Iterate time/FPS: %.3f ms/frame (%.1f FPS)", iterateTime * 1000, 1/iterateTime);

        ImGui::End();
    }
    // IMGui render
    ImGui::Render();
    SDL_SetRenderScale(sdlRenderer, io.DisplayFramebufferScale.x, io.DisplayFramebufferScale.y);
    SDL_RenderClear(sdlRenderer);

    // --- MY RENDERER ---
    double before = CycleTimer::currentSeconds();
    // Get the window surface
    windowSurface = SDL_GetWindowSurface(window);
    // Clear the surface
    SDL_ClearSurface(windowSurface, 0xFF, 0xFF, 0xFF, 0xFF);
    // Handle movement input
    eye.x += 0.1f * input.x;
    eye.z += 0.1f * input.y;
    cam.x += 0.1f * input.x;
    cam.z += 0.1f * input.y;
    renderer->SetCamera(eye, up, cam);
    // Render
    renderer->Render(windowSurface);
    // Flip the surface
    SDL_FlipSurface(windowSurface, SDL_FLIP_VERTICAL);
    // You gotta do this for some reason
    //SDL_UpdateWindowSurface(window);
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


void HandleInput(SDL_Event *event) {
    // if A or D is pressed
    // x = A + D
    // if W or S is pressed
    // y = W + S

    switch (event->key.key) {
        case SDLK_A:
            input.x = -1;
            break;
        case SDLK_D:
            input.x = 1;
            break;
        case SDLK_W:
            input.y = 1;
            break;
        case SDLK_S:
            input.y = -1;
            break;
        default:
            input.x = 0; input.y = 0;
            break;
    }
    // TODO: handle mouse
}

SDL_AppResult SDL_AppEvent(void *appstate, SDL_Event *event) {
    ImGui_ImplSDL3_ProcessEvent(event); // imgui event handler
    switch(event->type) {
        case SDL_EVENT_QUIT:
            return SDL_APP_FAILURE; // doing this will quit the program
            break;
        case SDL_EVENT_KEY_DOWN:
            HandleInput(event);
            SDL_Log("%d was pressed", event->key.scancode);
            break;
        default:
            break;
    }
    return SDL_APP_CONTINUE;
}

void SDL_AppQuit(void *appstate, SDL_AppResult result) {
    ImGui_ImplSDLRenderer3_Shutdown();
    ImGui_ImplSDL3_Shutdown();
    ImGui::DestroyContext();
    renderer->DestroyRenderer();
    SDL_DestroyRenderer(sdlRenderer);
    SDL_DestroyWindow(window);
    SDL_DestroySurface(windowSurface);
    SDL_Quit();
}
