//
// Created by iUV on 3/7/2025.
//
#include "TGAImage.h"
#include "geometry.h"

#ifndef SOFTWARERENDERER_MY_GL_HPP
#define SOFTWARERENDERER_MY_GL_HPP
void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color);
void line(Vec2i t0, Vec2i t1, TGAImage &image, TGAColor color);
void triangle(Vec3f *pts, TGAImage &image, float *zbuffer, TGAImage &texture, Vec2f texture_coords[3], int width);
void triangle(Vec3f *pts, TGAImage &image, float *zbuffer, TGAColor const &color, int width);
Matrix<float> homogonize(Vec3f v);
Vec3f dehomogonize(Matrix<float> const &m);
Vec3f project(Vec3f v, Matrix4x4f transfrom);
Vec3f rasterize(Vec3f v, Matrix4x4f m_viewport, Matrix4x4f m_proj, Matrix4x4f m_modelview);
Matrix4x4f LookAt(Vec3f eye, Vec3f center, Vec3f up);
Matrix4x4f world2screen(Vec3f v, int w, int h, float depth);
#endif //SOFTWARERENDERER_MY_GL_HPP
