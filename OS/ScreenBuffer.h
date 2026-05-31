//
// Created by iUV on 5/30/2026.
//

#ifndef SOFTWARERENDERER_SCREENBUFFER_H
#define SOFTWARERENDERER_SCREENBUFFER_H

class SDL_Texture;

/*
 * It will take in a SDL_Texture, lock it and manage the pixels exposed
 * by the SDL_LockTexture function for you to manage it in CPU
 * When this is deleted, it will relock the texture and destroy the texture object.
 * */

/// An abstraction for SDL_Texture to act like a framebuffer
/// Takes in a texture and exposes the Plot function
/// Texture and the pixel memory will be handled by this class

class ScreenBuffer {
private:
    // Internal texture
    SDL_Texture *mTexture;
    void *mPixel;
    int mWidth;
    int mHeight;
    int mPitch;
public:
    ScreenBuffer() = delete;
    ScreenBuffer(void** pixels, int width, int height, int pitch);
    ScreenBuffer(SDL_Texture *texture);
    ~ScreenBuffer();
    void Plot(int x, int y, int r, int g, int b, int a);
};


#endif //SOFTWARERENDERER_SCREENBUFFER_H
