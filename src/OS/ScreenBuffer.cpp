//
// Created by iUV on 5/30/2026.
//

#include "ScreenBuffer.h"
#include "SDL3/SDL.h"
#include <stdexcept>
#include <cstdint>

using namespace std;

ScreenBuffer::ScreenBuffer(SDL_Texture *texture): mTexture(texture)
{
    SDL_LockTexture(texture, nullptr, (void **) mPixel, &mPitch);
}

ScreenBuffer::~ScreenBuffer()
{
    SDL_UnlockTexture(mTexture);
}

void ScreenBuffer::Plot(int x, int y, int r, int g, int b, int a) {
    if(y * mPitch + x > mWidth * mHeight) {
        throw new std::out_of_range("Pixel out of range");
    }

    static_cast<uint32_t*>(mPixel)[y * mPitch + x] =
            a << 24 | r << 16 | g << 8 | b;
}