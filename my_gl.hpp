//
// Created by iUV on 3/7/2025.
//
#include "tgaimage.h"
#include "geometry.h"

#ifndef SOFTWARERENDERER_MY_GL_HPP
#define SOFTWARERENDERER_MY_GL_HPP

extern Matrix4x4f ModelView;
extern Matrix4x4f Projection;
extern Matrix4x4f Viewport;

/* Interface for both vertex and fragment shader*/
struct IShader {
    IShader() {
        uniform_M = Projection*ModelView;
        uniform_MIT = uniform_M;
        uniform_MIT.inverseTranspose();
    };
    Matrix4x4f uniform_M;
    Matrix<float> varying_uv = Matrix<float>(2, 3);
    Matrix4x4f uniform_MIT; // Invert transpose
    virtual ~IShader() = default;
    virtual Matrix<float> vertex(int iface, int nthvert) = 0;
    virtual bool fragment(Vec3f bar, TGAColor &color) = 0;
};

float *create_buffer(int width, int height);

void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color);
void line(Vec2i t0, Vec2i t1, TGAImage &image, TGAColor color);
void triangle(Vec3f *pts, TGAImage &image, float *zbuffer, int width, IShader &shader);
void wireframe_trig(Vec3f *pts, TGAImage &image, TGAColor color);
void triangle(Vec3f *pts, TGAImage &image, float *zbuffer, TGAColor const &color, int width);
Matrix<float> homogonize(Vec3f v, float h);
Vec3f dehomogonize(Matrix<float> const &m);
Vec3f project(Vec3f v);
Vec3f barycentric(Vec3f A, Vec3f B, Vec3f C, Vec3f P);
void Project(float coeff);
void LookAt(Vec3f eye, Vec3f center, Vec3f up);
void world2screen(Vec3f v, int w, int h, float depth);
void SetViewport(int width, int height, float depth);
void SetViewport(int x, int y, float w, float h, float depth);


#endif //SOFTWARERENDERER_MY_GL_HPP
