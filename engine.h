//
// Created by iUV on 5/12/2026.
//

#pragma once
#include "string"
#include "OS/Window.h"
constexpr float depth = 255.f;



class Engine {
private:
    void InitializeUniforms();
    Window window;
public:
    void Initialize();
    void Run();
    void RenderIMGUI();
    void HandleInput();
};


