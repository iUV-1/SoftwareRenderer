//
// Created by iuv on 7/30/25.
//


#ifndef RENDER_HPP
#define RENDER_HPP

//#include "Model.h
#include "tgaimage.h"
#include "geometry.h"

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
    //~Material();
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
    Scene(int width, int height, float depth, Vec3f &&light)
    : Width(width), Height(height), Depth(depth), Light(std::move(light)) {}
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
    Camera(const Vec3f &eye, const Vec3f &up, const Vec3f &cam): Eye(std::move(eye)), Up(std::move(up)), Cam(std::move(cam)) {}
    //~Camera();
    Vec3f Eye;
    Vec3f Up;
    Vec3f Cam;
};

class Renderer {
public:
    Renderer() { }
    //void LoadRendererSettings();
    void Render(SDL_Surface *surface);
    void DestroyRenderer();
    void SetCamera(const Vec3f &eye, const Vec3f &up, const Vec3f &cam);
    void SetupScene(int width, int height, float depth, Vec3f &&light);
    void SetupMaterial(std::string texPath, std::string normalPath, std::string specPath);
    void SetupModel(std::string modelPath);

    Model *mModel;
    Material *mMaterial;
    Scene *mScene;
    Camera *mCamera;

    Vec3f rasterize(IShader *shader, int iface, int nthvert);


};
#endif //RENDER_HPP
