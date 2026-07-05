//
// Created by iUV on 5/12/2026.
//

#pragma once

#include "../../OS/Window.h"
#include "../../Resources/Model.h"

class Scene;

// Interface for renderers
class IRenderer {
public:
    virtual ~IRenderer() = default;
    virtual void RenderIMGUI() = 0;
//    virtual void RenderModel(Model *model);
    virtual void RenderScene(Scene *scene) = 0;
};
