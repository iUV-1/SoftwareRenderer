//
// Created by iUV on 5/12/2026.
//

#pragma once

#include "../../OS/Window.h"
#include "../../Resources/Model.h"

// Interface for renderers
class IRenderer {
public:
    IRenderer() = delete;
    IRenderer(Window *window);
    ~IRenderer() = delete;
    virtual void RenderIMGUI();
    virtual void RenderModel(Model *model);
};