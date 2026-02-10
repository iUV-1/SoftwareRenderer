//
// Created by iuv on 7/30/25.
//


#ifndef RENDER_HPP
#define RENDER_HPP
#include <optional>
#include "tgaimage.h"
#include "geometry.h"
#include "RectBuffer.h"

class IShader;
class SDL_Surface;
class Model;

class Material {
public:
    //Material(const Material *mat);
    Material() { }
    Material(std::string texPath, std::string normalPath = "", std::string specPath = "")
    {
        TexFile->read_tga_file(texPath.c_str());
        TexFile->flip_vertically();
        if(normalPath != "") {
            NormalFile->read_tga_file(normalPath.c_str());
            NormalFile->flip_vertically();
            UseNormal = true;
        }
        if(specPath != "") {
            SpecularFile->read_tga_file(specPath.c_str());
            SpecularFile->flip_vertically();
            UseSpecular = true;
        }
    }
    ~Material() {
        delete NormalFile;
        delete TexFile;
        delete SpecularFile;
    }
    TGAImage *NormalFile = new TGAImage();
    TGAImage *TexFile = new TGAImage();
    TGAImage *SpecularFile = new TGAImage();
    bool UseNormal = false;
    bool UseSpecular = false;
};

class Scene {
public:
    Scene() { }
    Scene(const Scene *scene);
    Scene(int width, int height, float depth, const Vec3f &light)
    : Width(width), Height(height), Depth(depth), Light(light) {}
    //~Scene();
    int Width;
    int Height;
    float Depth;
    Vec3f Light;
};

class Camera {
public:
    Camera() { }
    // Copy constructor
    //Camera(const Camera *camera);
    Camera(const Vec3f &eye, const Vec3f &up, const Vec3f &cam, float fov): Eye(eye), Up(up), Cam(cam), Fov(fov) {}
    //~Camera();
    Vec3f Eye;
    Vec3f Up;
    Vec3f Cam;
    float Fov;
};

class Renderer {
public:
    Renderer() { }
    ~Renderer();
    //void LoadRendererSettings();
    void Render(SDL_Surface *surface);
    void DestroyRenderer();
    void SetCamera(const Vec3f &eye, const Vec3f &up, const Vec3f &cam, float fov);
    void SetupScene(int width, int height, float depth, const Vec3f &light);
    void SetupMaterial(std::string texPath, std::string normalPath, std::string specPath);
    void SetupModel(std::string modelPath);
    void ChangeResolution(int width, int height);
    void SetupBuffers();

    Model *mModel;
    Material *mMaterial;
    Scene *mScene;
    Camera *mCamera;
    RectBuffer *mShadowBuffer;
    RectBuffer *mZbuffer;
private:
    Vec3f rasterize(IShader *shader, int iface, int nthvert);
};
#endif //RENDER_HPP
