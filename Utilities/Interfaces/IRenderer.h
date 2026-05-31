//
// Created by iUV on 5/12/2026.
//

#ifndef SOFTWARERENDERER_IRENDERER_H
#define SOFTWARERENDERER_IRENDERER_H

// Interface for renderers
class IRenderer {
public:
    virtual void Initialize();
    virtual void RenderIMGUI();
    virtual void RenderMesh();
    virtual void End();
};


#endif //SOFTWARERENDERER_IRENDERER_H
