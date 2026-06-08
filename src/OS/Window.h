//
// Created by iUV on 5/27/2026.
//

#ifndef SOFTWARERENDERER_WINDOW_H
#define SOFTWARERENDERER_WINDOW_H
#include "imgui.h"
#include "SDL3/SDL.h"
#include <string>
#include "ScreenBuffer.h"
#include "../../Utilities/Math/geometry.h"
/// Handler for SDL3 and initialize ImGUI
class Window{
private:

    std::string mTitle;
    SDL_Window *mWindow;
    SDL_Renderer *mSDLRenderer;
    SDL_Surface *mSurface;
//    /// What will be rendered
//    ScreenBuffer *frontBuffer;
//    ScreenBuffer *backBuffer;
    /// Double buffer
    //SDL_Texture *backTexture;
    ImGuiIO io;
    Vec2f mRelativeMouse;
    const bool *mKeyboardState;
public:
    int mWidth;
    int mHeight;
    Window(std::string title, int width, int height);
    ~Window();
    void Plot(int x, int y, int r, int g, int b, int a);
    void Update();
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
};


#endif //SOFTWARERENDERER_WINDOW_H
