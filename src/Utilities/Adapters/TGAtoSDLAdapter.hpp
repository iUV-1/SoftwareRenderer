//
// Created by iUV on 1/6/2026.
//


#ifndef SOFTWARERENDERER_TGATOSDLADAPTER_HPP
#define SOFTWARERENDERER_TGATOSDLADAPTER_HPP

#include "../../Resources/tgaimage.h"
#include "SDL3/SDL.h"
/// Adapter for TGA to SDL
namespace TGAtoSDLAdapter {
    void static SetTGAPixel(SDL_Surface *surface, int x, int y, TGAColor color) {
        if(surface->format != SDL_PIXELFORMAT_XRGB8888) {
            // bruh
            throw std::invalid_argument("Unsupported format");
        }
        Uint32* target_pixel = (Uint32*)surface->pixels + (y * surface->pitch/4 + x);

        int a = color.a;
        int r = color.r;
        int g = color.g;
        int b = color.b;
        *target_pixel = a << 24 | r << 16 | g << 8 | b;
    }
}


#endif //SOFTWARERENDERER_TGATOSDLADAPTER_HPP
