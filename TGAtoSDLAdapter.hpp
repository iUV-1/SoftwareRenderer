//
// Created by iUV on 1/6/2026.
//


#ifndef SOFTWARERENDERER_TGATOSDLADAPTER_HPP
#define SOFTWARERENDERER_TGATOSDLADAPTER_HPP

#include "tgaimage.h"
#include "SDL3/SDL.h"
/// Adapter for TGA to SDL
///
class TGAtoSDLAdapter {
public:
    void static SetTGAPixel(SDL_Surface *surface, int x, int y, TGAColor color);
};


#endif //SOFTWARERENDERER_TGATOSDLADAPTER_HPP
