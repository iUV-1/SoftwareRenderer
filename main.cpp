//
// Created by iUV on 1/6/2026.
//
#define SDL_MAIN_USE_CALLBACKS 1

#include "pch.h"
#include <SDL3/SDL_main.h>
#include <string>

#include "render.hpp"
#include "CycleTimer.hpp"

SDL_Window *window;
SDL_Surface *windowSurface;
SDL_Renderer *sdlRenderer;

constexpr int width = 640;
constexpr int height = 640;
constexpr float depth = 255.f;
Renderer *renderer;

Vec2i input;
Vec3f eye(1, 1, 3);
Vec3f up(0, 1, 0);
Vec3f cam(0, 0, 0);

SDL_AppResult SDL_AppInit(void **appstate, int argc, char *argv[]) {
    if(!SDL_Init(SDL_INIT_VIDEO | SDL_INIT_JOYSTICK)) {
        return SDL_APP_FAILURE;
    }
    window = SDL_CreateWindow("Software Renderer", width, height, SDL_WINDOW_RESIZABLE);
    // Initialize the renderer
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

    sdlRenderer = SDL_CreateRenderer(window, nullptr);
    return SDL_APP_CONTINUE;
}

SDL_AppResult SDL_AppIterate(void *appstate) {
    SDL_RenderClear(sdlRenderer);
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
    // Blit
    //SDL_BlitSurface(bmpSurface, nullptr, windowSurface, nullptr);
    // You gotta do this for some reason
    SDL_UpdateWindowSurface(window);
    double timeTaken = CycleTimer::currentSeconds() - before;
    double fps = 1 / timeTaken;
    std::string fpsStr = "FPS: " + std::to_string(fps);
    SDL_Log("%.2f ms - %.1f FPS\n", 1000.f * timeTaken, fps);

    SDL_RenderDebugText(sdlRenderer, 10, 10, fpsStr.c_str());
    //SDL_RenderPresent(sdlRenderer);
    return SDL_APP_CONTINUE;
}


void HandleInput(SDL_Event *event) {
    // if A or D is pressed
    // x = A + D
    // if W or S is pressed
    // y = W + S

    switch (event->key.key) {
        case SDLK_A:
            input.x -= 1;
            break;
        case SDLK_D:
            input.x += 1;
            break;
        case SDLK_W:
            input.y += 1;
            break;
        case SDLK_S:
            input.y -= 1;
            break;
        default:
            input.x = 0; input.y = 0;
            break;
    }
    // TODO: handle mouse
}

SDL_AppResult SDL_AppEvent(void *appstate, SDL_Event *event) {
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
    renderer->DestroyRenderer();
    SDL_Quit();
}
