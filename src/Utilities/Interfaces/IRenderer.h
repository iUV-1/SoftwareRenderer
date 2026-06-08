//
// Created by iUV on 5/12/2026.
//

#pragma once

#include "../../OS/Window.h"
#include "../../Resources/Model.h"

// Interface for renderers
class IRenderer {
public:
    IRenderer(Window *window);

    virtual void RenderIMGUI();
//    virtual void RenderModel(Model *model);
    virtual void RenderScene(Scene *scene);
};
