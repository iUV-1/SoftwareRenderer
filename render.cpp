//
// Created by iuv on 7/30/25.
//
#include <string>
#include "pch.h"

#include "render.hpp"
#include "model.h"
#include "tgaimage.h"
#include "my_gl.hpp"
//#include "shaders.hpp"

#include "shaders/PhongShadowShader.h"
#include "shaders/DepthShaderImage.h"

using namespace std;

Vec3f Renderer::rasterize(IShader *shader, int iface, int nthvert) {
    // Apply vertex shader
    Vec4f homogonized = shader->vertex(iface, nthvert);
    Vec3f result = dehomogonize(homogonized);
    // Round the result to apply to screen
    result.x = std::round(result.x);
    result.y = std::round(result.y);
    result.z = std::round(result.z);

    return result;
}

void Renderer::SetupModel(string modelPath) {
    mModel = new Model(modelPath.c_str());
}

void Renderer::SetupMaterial(std::string texPath, std::string normalPath, std::string specPath) {
    mMaterial = new Material(texPath, normalPath, specPath);
}

void Renderer::SetupScene(int width, int height, float depth, const Vec3f &light) {
    mScene = new Scene(width, height, depth, light);
}

void Renderer::SetCamera(const Vec3f &eye, const Vec3f &up, const Vec3f &cam) {
    mCamera = new Camera(eye, up, cam);
}

void Renderer::DestroyRenderer() {
    delete mModel;
    delete mMaterial;
    delete mScene;
    delete mCamera;
}

void Renderer::Render(SDL_Surface *surface) {
    /* Depth map */
    // Init depth buffer
    //float *depth_buffer_arr = create_buffer(mScene->Width, mScene->Height);
    RectBuffer depthBufferArr(mScene->Width, mScene->Height);
    // Init shader
    LookAt(mScene->Light, mCamera->Cam, mCamera->Up); // Render from the mLight
    Project(0); // Render mLight in orthographic mode
    SetViewport(mScene->Width / 8, mScene->Height/8, mScene->Width * 3./4, mScene->Height * 3./4, mScene->Depth); // Clamp the image into the center with margins (3/4 of the screen)

    auto depth_shader = DepthShaderImage(mModel, mScene);
    // Render
    for(int i = 0; i < mModel->nfaces(); ++i) {
        Vec3f screen_coords[3];
        for (int j=0; j<3; ++j)
            screen_coords[j] = rasterize(&depth_shader, i, j);

        Vec3f n = (screen_coords[2]-screen_coords[0])^(screen_coords[1]-screen_coords[0]);
        n.normalize();
        float view_dir_intensity = mCamera->Eye*n;

        if (view_dir_intensity<1)
            // Doesn't save the image since we don't need it
            triangle(screen_coords, mScene->Width, mScene->Height, depthBufferArr);
    }
    Matrix4x4f M_Shadow = Viewport*Projection*ModelView;

    /* SSAO */
    /// Get depth from Camera
    //auto depth_from_cam = TGAImage(mScene->Width, mScene->Height, TGAImage::RGB);
    // Setup GL
    LookAt(mCamera->Eye, mCamera->Cam, mCamera->Up);
    Project(-1/(mCamera->Eye-mCamera->Cam).norm());
    //SetViewport(mWidth/8, mHeight/8, mWidth*3./4, mHeight*3./4, mDepth);

    // Setup zbuffer
/*    float *depth_cam_buf = create_buffer(mWidth, mHeight);
    auto depth_shader_image = DepthShaderImage();
    for(int i = 0; i < model->nfaces(); ++i) {
        Vec3f screen_coords[3];
        for (int j=0; j<3; ++j)
            screen_coords[j] = rasterize(&depth_shader_image, i, j);

        Vec3f n = (screen_coords[2]-screen_coords[0])^(screen_coords[1]-screen_coords[0]);
        n.normalize();
        float view_dir_intensity = mEye*n;

        if (view_dir_intensity<1)
            triangle(screen_coords, depth_from_cam, depth_cam_buf, mWidth, depth_shader_image);
    }
    depth_from_cam.flip_vertically();*/

    /// SSAO
    /*auto ssao_frame = TGAImage(mWidth, mHeight, TGAImage::RGB);

    for (int x=0; x<mWidth; x++) {
        for (int y=0; y<mHeight; y++) {
            if (depth_cam_buf[x+y*mWidth] < -1e5) continue;
            float total = 0;
            for (float a=0; a<M_PI*2-1e-4; a += M_PI/4) {
                total += M_PI/2 - max_elevation_angle(depth_cam_buf, Vec2f(x, y), Vec2f(cos(a), sin(a)));
            }
            total /= (M_PI/2)*8;
            total = pow(total, 100.f);
            ssao_frame.set(x, y, TGAColor(total*255, total*255, total*255));
        }
    }
    ssao_frame.flip_vertically();*/
    /* Render */
    auto frame = TGAImage(mScene->Width, mScene->Height, TGAImage::RGB);
    // Setup zbuffer
//    float *zbuffer = create_buffer(mScene->Width, mScene->Height);
    RectBuffer zbuffer(mScene->Width, mScene->Height);
    // GouraudShader shader = GouraudShader();
    // PhongShader shader = PhongShader();
    // GouraudShaderReference shader = GouraudShaderReference();
    Matrix4x4f MVP = Viewport*Projection*ModelView;
    MVP.invert();
    PhongShaderShadow shader = PhongShaderShadow(mMaterial, mCamera, mModel, mScene, M_Shadow, depthBufferArr);
    //PhongShader shader = PhongShader();

    for (int i=0; i< mModel->nfaces(); ++i) {
        Vec3f screen_coords[3];

        for (int j=0; j<3; ++j)
            screen_coords[j] = rasterize(&shader, i, j);

        // calculate normal
        // ^ is an overloaded operator that performs cross product calculation
        // world_coords[2] - world_coords[0] and the other are 2 vectors pointing from point
        // world_coords[0].
        Vec3f n = (screen_coords[2]-screen_coords[0])^(screen_coords[1]-screen_coords[0]);
        n.normalize();
        // calculate mEye intensity by dot product between normal and mEye vector
        float view_dir_intensity = mCamera->Eye*n;
        // back face culling

        if (view_dir_intensity<1) {
            //triangle(screen_coords, frame, zbuffer, mScene->Width, shader);
            triangle(screen_coords, surface, zbuffer, shader);
            //wireframe_trig(screen_coords, frame, TGAColor(255, 255, 255, 255));
        }
    }
    // set origin to the bottom left corner
    // frame.flip_vertically();
}
