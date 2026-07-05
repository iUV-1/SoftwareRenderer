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

#include "Renderer/SoftwareRenderer/render.hpp"
#include "Utilities/CycleTimer.hpp"
#include "engine.h"


int width = 800;
int height = 600;
constexpr float depth = 255.f;
Engine *engine;
struct GlobalAppState {
    const bool *keystate;
    ImGuiIO imGuiIo;
};

SDL_AppResult SDL_AppInit(void **appstate, int argc, char *argv[]) {
    // Get the global struct by casting
    engine = new Engine();
    engine->Initialize(argc, argv);
    return SDL_APP_CONTINUE;
}

SDL_AppResult SDL_AppIterate(void *appstate) {
    engine->Update();
    return SDL_APP_CONTINUE;
}

SDL_AppResult SDL_AppEvent(void *appstate, SDL_Event *event) {
    ImGui_ImplSDL3_ProcessEvent(event); // imgui event handler

    switch(event->type) {
        case SDL_EVENT_QUIT:
            return SDL_APP_FAILURE; // doing this will quit the program
            break;
        case SDL_EVENT_KEY_DOWN:
            break;
        case SDL_EVENT_WINDOW_RESIZED:
            engine->ResizeWindow(event->window.data1, event->window.data2);
            break;
        default:
            break;
    }
    return SDL_APP_CONTINUE;
}

void SDL_AppQuit(void *appstate, SDL_AppResult result) {
    delete (GlobalAppState*)appstate;
    delete engine;
    //SDL_DestroySurface(windowSurface);
    SDL_Quit();
}
