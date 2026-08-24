//
// Created by iuv on 7/30/25.
//
#include <string>
#include "../../../pch.h"

#include "render.hpp"
#include "../../Resources/Mesh.h"
#include "../../Resources/Scene.h"
#include "../../Resources/tgaimage.h"
#include "../../Utilities/Math/math.h"
#include "my_gl.hpp"

#include "./shaders/PhongShadowShader.h"
#include "./shaders/DepthShaderImage.h"

using namespace std;

Vec3f Renderer::rasterize(IShader *shader, int iface, int nthvert)
{
    // Apply vertex shader
    Vec4f homogonized = shader->vertex(iface, nthvert);
    Vec3f result = dehomogonize(homogonized);
    // Round the result to apply to screen
//    result.x = round(result.x);
//    result.y = round(result.y);
//    result.z = round(result.z);

    return result;
}

void Renderer::SetupBuffers(int width, int height)
{
    // Invalidate the previous buffer (if it exists)
    // Used to change the resolution of the screen and camera
    if(mShadowBuffer) delete mShadowBuffer;
    if(mZbuffer) delete mZbuffer;
    mShadowBuffer = new RectBuffer(width, height);
    mZbuffer = new RectBuffer(width, height);
}

void Renderer::SetupUniforms(Matrix4x4f const &projection, Matrix4x4f const &modelview, Matrix4x4f const &viewport)
{
    mProjection = projection;
    mModelView = modelview;
    mViewport = viewport;
}

Renderer::~Renderer()
{
    delete mCamera;
    delete mShadowBuffer;
    delete mZbuffer;
}

void Renderer::RenderIMGUI()
{

}

void Renderer::RenderScene(Scene *scene)
{
    /* Depth map */
    mScene = scene;
    mCamera = scene->GetCamera();
    Model *model = mScene->GetModels()[0];
    Vec3f light = *mScene->GetLight();
    // Init depth buffer
    mShadowBuffer->resetBuffer();
    // Init shader
    Matrix4x4f shadowModelView = LookAt(light, mCamera->Cam, mCamera->Up);
    Matrix4x4f shadowProjection = OrthoProject(0); // Render mLight in orthographic mode
    auto depth_shader = DepthShaderImage(mViewport, shadowProjection, shadowModelView, model, mCamera->Depth);
    // Render
    for(int i = 0; i < model->mMesh.nfaces(); ++i) {
        Vec3f screen_coords[3];
        for (int j=0; j<3; ++j)
            screen_coords[j] = rasterize(&depth_shader, i, j);
        // calculate normal
        // ^ is an overloaded operator that performs cross product calculation
        // world_coords[2] - world_coords[0] and the other are 2 vectors pointing from point
        // world_coords[0].
        Vec3f n = (screen_coords[2]-screen_coords[0])^(screen_coords[1]-screen_coords[0]);
        n.normalize();

        float view_dir_intensity = mCamera->Eye*n;

        if (view_dir_intensity<1)
            // Doesn't save the image since we don't need it
            triangle(screen_coords, mCamera->Width, mCamera->Height, *mShadowBuffer);
    }
    Matrix4x4f M_Shadow = mViewport*shadowProjection*mModelView;

    /* SSAO */
    /// Get depth from Camera
    //auto depth_from_cam = TGAImage(mScene->Width, mScene->Height, TGAImage::RGB);
    // Setup GL
//    LookAt(mCamera->Eye, mCamera->Cam, mCamera->Up);
//    Project(mCamera->Fov, mScene->Width, mScene->Height, 0.1, mScene->Depth);

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
    mZbuffer->resetBuffer();
    // GouraudShader shader = GouraudShader();
    // PhongShader shader = PhongShader();
    // GouraudShaderReference shader = GouraudShaderReference();
//    Matrix4x4f MVP = mViewport*mProjection*mModelView;
//    Matrix4x4f MVP = mModelView * mViewport * mProjection;
//    MVP.invert();
    PhongShaderShadow shader = PhongShaderShadow(mViewport, mProjection, mModelView, model, light, M_Shadow, *mShadowBuffer);
    //PhongShader shader = PhongShader();
    for (int i=0; i< model->mMesh.nfaces(); ++i) {
        Vec3f screen_coords[3];

        for (int j=0; j<3; ++j)
            screen_coords[j] = rasterize(&shader, i, j);


        Vec3f n = (screen_coords[2]-screen_coords[0])^(screen_coords[1]-screen_coords[0]);
        n.normalize();
        // calculate mEye intensity by dot product between normal and mEye vector
        float view_dir_intensity = mCamera->Eye*n;
        // back face culling
        if (view_dir_intensity<1) {
            //triangle(screen_coords, frame, zbuffer, mScene->Width, shader);
            triangle(screen_coords, mWindow, *mZbuffer, shader);
            //wireframe_trig(screen_coords, frame, TGAColor(255, 255, 255, 255));
        }
    }
}