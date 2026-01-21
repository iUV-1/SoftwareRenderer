//
// Created by iUV on 3/7/2025.
//


#pragma once
#include "tgaimage.h"
#include "geometry.h"
#include "model.h"
#include "render.hpp"

extern Matrix4x4f ModelView;
extern Matrix4x4f Projection;
extern Matrix4x4f Viewport;

#define TEX2D(tex, uv) (tex->get(uv.u * tex->get_width(), uv.v * tex->get_height()))

class RectBuffer {
public:
    RectBuffer() = delete;
    RectBuffer(int const w, const int h): width(w), size(w*h) {
        data = new float[size];
    }
    ~RectBuffer() {
        delete[] data;
    }

    float &operator[] (size_t const i) {
        if (i >= size) throw std::out_of_range("Buffer index out of range");
        return data[i];
    }
    size_t width;
    size_t size;
    float *data;
};

/* Interface for both vertex and fragment shader*/
class IShader {
public:
    IShader() {}
    IShader(Scene *scene) {
        uniform_M = Projection*ModelView;
        uniform_MIT = uniform_M;
        uniform_MIT.inverseTranspose();
        width = scene->Width;
    }
    IShader(Material *mat, Camera *cam, Model *model, Scene *scene)
            : uniform_Material(mat), uniform_Camera(cam), uniform_Model(model), uniform_Scene(scene){
//    :width(w){
        //input = i;
        uniform_M = Projection*ModelView;
        uniform_MIT = uniform_M;
        uniform_MIT.inverseTranspose();
        width = scene->Width;
    };
    //IShader() = default; // don't use this

    Matrix4x4f uniform_M;
    Matrix<float> varying_uv = Matrix<float>(2, 3);
    Matrix4x4f uniform_MIT; // Invert transpose
    Material *uniform_Material = nullptr;
    Camera  *uniform_Camera = nullptr;
    Model *uniform_Model = nullptr;
    Scene *uniform_Scene = nullptr;
    int width;
    //virtual ~IShader() = default;
    virtual Vec4f vertex(int iface, int nthvert) = 0;
    virtual bool fragment(Vec3f bar, TGAColor &color) = 0;
};

float *create_buffer(int width, int height);

void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color);
void line(Vec2i t0, Vec2i t1, TGAImage &image, TGAColor color);
void triangle(Vec3f *pts, TGAImage &image, float *zbuffer, int width, IShader &shader);
void wireframe_trig(Vec3f *pts, TGAImage &image, TGAColor color);
void triangle(Vec3f *pts, TGAImage &image, float *zbuffer, TGAColor const &color, int width);
void triangle(Vec3f *pts, SDL_Surface *surface, float *zbuffer, int width, IShader &shader);
void triangle(Vec3f *pts, SDL_Surface *surface, RectBuffer &zbuffer, IShader &shader);
void triangle(Vec3f *pts, int w, int h, RectBuffer &zbuffer);
//Matrix<float> homogonize(Vec3f v, float h);
//Vec3f dehomogonize(Matrix<float> const &m);
Vec3f project(Vec3f v);
Vec3f barycentric(Vec3f A, Vec3f B, Vec3f C, Vec3f P);
void Project(float coeff);
void LookAt(Vec3f eye, Vec3f center, Vec3f up);
void world2screen(Vec3f v, int w, int h, float depth);
void SetViewport(int width, int height, float depth);
void SetViewport(int x, int y, float w, float h, float depth);

