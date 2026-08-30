//
// Created by iUV on 5/27/2026.
//

#ifndef SOFTWARERENDERER_WINDOW_H
#define SOFTWARERENDERER_WINDOW_H
#include "imgui.h"
#include "SDL3/SDL.h"
#include <string>
#include "../../Utilities/Math/geometry.h"
/// Handler for SDL3 and initialize ImGUI
class Window{
private:
    std::string mTitle;
    SDL_Window *mWindow;
    SDL_Renderer *mSDLRenderer;
    SDL_Surface *mSurface;
    SDL_Texture *mTex;
//    /// What will be rendered
//    ScreenBuffer *frontBuffer;
//    ScreenBuffer *backBuffer;
    /// Double buffer
    //SDL_Texture *backTexture;
    ImGuiIO io;
    Vec2f mRelativeMouse;
    const bool *mKeyboardState;
    bool mMouseState = true;
public:
    int mWidth;
    int mHeight;
    Window(std::string title, int width, int height);
    ~Window();
    void Plot(int x, int y, int r, int g, int b, int a);
    void Update();
    void StartUpdate();
    void EndUpdate();
    void Resize(int width, int height);
    /// The global state of the input.
    /// Technically not handled by the window but since
    /// the backend is in SDL3, the window hold it for clarity
    const bool* GetKeyboardState()
    {
        return mKeyboardState;
    };
    const Vec2f* GetRelativeMouse()
    {
        return &mRelativeMouse;
    }

    void ToggleRelativeMouse();
};


#endif //SOFTWARERENDERER_WINDOW_H
