//
// Created by iUV on 3/7/2025.
//


#pragma once
#include "../../Resources/tgaimage.h"
#include "../../Utilities/Math/geometry.h"
#include "../../Resources/model.h"
#include "render.hpp"

extern Matrix4x4f ModelView;
extern Matrix4x4f Projection;
extern Matrix4x4f Viewport;

#define TEX2D(tex, uv) (tex->get(uv.u * tex->get_width(), uv.v * tex->get_height()))

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
void Project(float fov, float width, float height, float near, float far);
void LookAt(Vec3f eye, Vec3f center, Vec3f up);
void world2screen(Vec3f v, int w, int h, float depth);
void SetViewport(int width, int height, float depth);
void SetViewport(int x, int y, float w, float h, float depth);