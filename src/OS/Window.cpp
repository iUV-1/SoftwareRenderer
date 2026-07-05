//
// Created by iUV on 5/27/2026.
//

#include "Window.h"
#include "SDL3/SDL.h"
#include "backends/imgui_impl_sdlrenderer3.h"
#include "backends/imgui_impl_sdl3.h"
#include <stdexcept>

Window::Window(std::string title, int width, int height): mTitle(title), mWidth(width), mHeight(height) {
    mWindow = SDL_CreateWindow(title.c_str(), width, height, SDL_WINDOW_RESIZABLE);
    SDL_SetWindowRelativeMouseMode(mWindow, true);

    mSDLRenderer = SDL_CreateRenderer(mWindow, nullptr);
    if (mSDLRenderer == nullptr)
    {
        SDL_Log("Error: SDL_CreateRenderer(): %s\n", SDL_GetError());
        throw std::invalid_argument("Failed to create SDL renderer!");
    }

    mKeyboardState = SDL_GetKeyboardState(nullptr);
//    frontBuffer = new ScreenBuffer(SDL_CreateTextureFromSurface(sdlRenderer, surf));
//    backBuffer = new ScreenBuffer(SDL_CreateTextureFromSurface(sdlRenderer, surf));
    // --- SETUP IMGUI ---
    // imgui stuff
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    io = ImGui::GetIO();    (void)io;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
    // Setup Dear ImGui style
    ImGui::StyleColorsDark();
    //ImGui::StyleColorsLight();
    float main_scale = SDL_GetDisplayContentScale(SDL_GetPrimaryDisplay());

    // Setup scaling
    ImGuiStyle& style = ImGui::GetStyle();
    style.ScaleAllSizes(main_scale);        // Bake a fixed style scale. (until we have a solution for dynamic style scaling, changing this requires resetting Style + calling this again)
    style.FontScaleDpi = main_scale;        // Set initial font scale. (using io.ConfigDpiScaleFonts=true makes this unnecessary. We leave both here for documentation purpose)
    // Setup Platform/Renderer backends
    ImGui_ImplSDL3_InitForSDLRenderer(mWindow, mSDLRenderer);
    ImGui_ImplSDLRenderer3_Init(mSDLRenderer);
}

/// No guard checking for performance. Good luck!
void Window::Plot(int x, int y, int r, int g, int b, int a) {
    static_cast<uint32_t*>(mSurface->pixels)[y * mWidth + x] =
            a << 24 | r << 16 | g << 8 | b;
}

void Window::StartUpdate()
{
    mSurface = SDL_GetWindowSurface(mWindow);

    ImGui_ImplSDLRenderer3_NewFrame();
    ImGui_ImplSDL3_NewFrame();
    ImGui::NewFrame();

    // Input
    SDL_GetRelativeMouseState(&mRelativeMouse.x, &mRelativeMouse.y);
    // Swap surfaces
    SDL_RenderClear(mSDLRenderer);
    SDL_ClearSurface(mSurface, 0xff, 0xff, 0xff, 0xff);
}

void Window::EndUpdate()
{
    SDL_UpdateWindowSurface(mWindow);

    // ImGUI Render
    ImGui::Render();
    SDL_SetRenderScale(mSDLRenderer, io.DisplayFramebufferScale.x, io.DisplayFramebufferScale.y);
    ImGui_ImplSDLRenderer3_RenderDrawData(ImGui::GetDrawData(), mSDLRenderer);

    SDL_RenderPresent(mSDLRenderer);
}

void Window::Update() {



}

Window::~Window() {
    ImGui_ImplSDLRenderer3_Shutdown();
    ImGui_ImplSDL3_Shutdown();
    ImGui::DestroyContext();
    SDL_DestroyRenderer(mSDLRenderer);
    SDL_DestroyWindow(mWindow);
    //SDL_DestroySurface(windowSurface);
    SDL_Quit();
}