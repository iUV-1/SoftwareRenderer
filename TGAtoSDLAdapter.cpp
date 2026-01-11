//
// Created by iUV on 1/6/2026.
//

#include "TGAtoSDLAdapter.hpp"

void TGAtoSDLAdapter::SetTGAPixel(SDL_Surface *surface, int x, int y, TGAColor color) {
    if(surface->format != SDL_PIXELFORMAT_RGBA32) {
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