#include "tgaimage.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <limits>
#include "model.h"
#include "geometry.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
const TGAColor green = TGAColor(0,   255, 0,  255);

#define Vec2i Vec2<int>

void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color);
void line(Vec2i t0, Vec2i t1, TGAImage &image, TGAColor color);
void triangle(Vec3f *pts, TGAImage &image, float *zbuffer, TGAImage &texture, Vec2f texture_coords[3]);
void triangle(Vec3f *pts, TGAImage &image, float *zbuffer, TGAColor const &color);

Vec3f world2screen(Vec3f v);
Model *model = nullptr;
const int width = 800;
const int height = 800;

Vec3f world2screen(Vec3f v) {
    return Vec3f(static_cast<int>((v.x+1.)*width/2.+.5), static_cast<int>((v.y+1.)*height/2. + .5), v.z);
}

int main(int argc, char** argv) {
    if(argc < 2) {
        model = new Model("obj/african_head.obj");
    } else {
        model = new Model(argv[1]);
    }

    TGAImage tex_file(1024,1024,TGAImage::RGB);
    if(argc < 3) {
        tex_file.read_tga_file("obj/UV Grid.tga");
    } else {
        tex_file.read_tga_file(argv[2]);
    }

    TGAImage frame(width, height, TGAImage::RGB);
    auto *zbuffer = new float[width * height];
    for (int i=width*height; i--; zbuffer[i] = -std::numeric_limits<float>::max());
    Vec3f light(0,0, -1);
    for (int i=0; i<model->nfaces(); ++i) {
        std::vector<int> face = model->face(i);
        std::vector<int> face_tex = model->face_tex(i);
        Vec3f screen_coords[3];
        Vec3f world_coords[3];
        Vec2f texture_coords[3];
        for (int j=0; j<3; ++j) {
            Vec3f v = model->vert(face[j]);
            texture_coords[j] = model->texcoord(face_tex[j]);
            screen_coords[j] = world2screen(v);
            world_coords[j] = v;
        }
        // calculate normal
        // ^ is an overloaded operator that performs cross product calculation
        // world_coords[2] - world_coords[0] and the other are 2 vectors pointing from point
        // world_coords[0].
        Vec3f n = (world_coords[2]-world_coords[0])^(world_coords[1]-world_coords[0]);
        n.normalize();
        // calculate light intensity by dot product between normal and light vector
        float intensity = n*light;
        float view_dir_intensity = n*Vec3f(0, 0, -1);
        // back face culling
        triangle(screen_coords, frame, zbuffer, tex_file, texture_coords);
        //triangle(screen_coords, frame, zbuffer, TGAColor(intensity*255, intensity*255, intensity*255, 255));

        if (view_dir_intensity>0) {
        }
    }
    frame.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    frame.write_tga_file("./output.tga");
    delete[] zbuffer;
    delete model;
    return 0;
}


void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color) {
    bool steep = false;
    // If it's too steep, swap it due to the algo struggling with steep lines
    if(std::abs(x0-x1) < std::abs(y0-y1)) {
        std::swap(x0, y0);
        std::swap(x1,y1);
        steep = true;
    }

    if(x0>x1) {
        std::swap(x0, x1);
        std::swap(y0, y1);
    }

    int dx = x1 - x0;
    int dy = y1 - y0;
    // slope of the line
    float derror = std::abs(dy)*2;
    float error = 0;
    int y = y0;

    for (float x=x0; x<=x1; x++) {
        if(steep) {
            image.set(y, x, color);
        } else {
            image.set(x, y, color);
        }
        // if the error gets too big, increment y
        error += derror;
        if (error>dx) {
            y += (y1>y0?1:-1);
            error -= dx*2;
        }
    }
}

// Overload with vec2
void line(Vec2<int> t0, Vec2<int> t1, TGAImage &image, TGAColor color) {
    int x0 = t0.x;
    int y0 = t0.y;
    int x1 = t1.x;
    int y1 = t1.y;
    bool steep = false;
    // If it's too steep, swap it due to the algo struggling with steep lines
    if(std::abs(x0-x1) < std::abs(y0-y1)) {
        std::swap(x0, y0);
        std::swap(x1,y1);
        steep = true;
    }

    if(x0>x1) {
        std::swap(x0, x1);
        std::swap(y0, y1);
    }

    int dx = x1 - x0;
    int dy = y1 - y0;
    // slope of the line
    float derror = std::abs(dy)*2;
    float error = 0;
    int y = y0;

    for (float x=x0; x<=x1; x++) {
        if(steep) {
            image.set(y, x, color);
        } else {
            image.set(x, y, color);
        }
        // if the error gets too big, increment y
        error += derror;
        if (error>dx) {
            y += (y1>y0?1:-1);
            error -= dx*2;
        }
    }
}

Vec3f barycentric(Vec3f A, Vec3f B, Vec3f C, Vec3f P) {
    Vec3f s[2];
    for (int i=2; i--; ) {
        s[i][0] = C[i]-A[i];
        s[i][1] = B[i]-A[i];
        s[i][2] = A[i]-P[i];
    }
    Vec3f u = s[0]^s[1];

    if (std::abs(u[2])>1e-2) // dont forget that u[2] is integer. If it is zero then triangle ABC is degenerate

        return Vec3f(1.f-(u.x+u.y)/u.z, u.y/u.z, u.x/u.z);

    return Vec3f(-1,1,1); // in this case generate negative coordinates, it will be thrown away by the rasterizator}
}

void triangle(Vec3f *pts, TGAImage &image, float *zbuffer, TGAImage &texture, Vec2f texture_coords[3]) {
    Vec2f bboxmin( std::numeric_limits<float>::max(),  std::numeric_limits<float>::max());
    Vec2f bboxmax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    Vec2f clamp(image.get_width()-1, image.get_height()-1);
    for (int i=0; i<3; i++) {
        for (int j=0; j<2; j++) {
            bboxmin[j] = std::max(0.f,      std::min(bboxmin[j], pts[i][j]));
            bboxmax[j] = std::min(clamp[j], std::max(bboxmax[j], pts[i][j]));
        }
    }
    Vec3f P;
    for (P.x=bboxmin.x; P.x<=bboxmax.x; ++P.x) {
        for (P.y=bboxmin.y; P.y<=bboxmax.y; ++P.y) {
            Vec3f bc_screen  = barycentric(pts[0], pts[1], pts[2], P);
            if (bc_screen.x<0 || bc_screen.y<0 || bc_screen.z<0) continue;
            P.z = 0;
            for (int i=0; i<3; ++i) {
                P.z += pts[i][2]*bc_screen[i];
            }
            auto idx = static_cast<size_t>(P.x + P.y * width);
            float u = bc_screen.x * texture_coords[0].u + bc_screen.y * texture_coords[1].u + bc_screen.z * texture_coords[2].u;
            float v = bc_screen.x * texture_coords[0].v + bc_screen.y * texture_coords[1].v + bc_screen.z * texture_coords[2].v;
            if(zbuffer[idx] < P.z) {
                zbuffer[idx] = P.z;
                image.set(P.x, P.y, texture.get(u*texture.get_width(), v*texture.get_height()));
            }
        }
    }
}
// Overload that is just intensity
void triangle(Vec3f *pts, TGAImage &image, float *zbuffer, TGAColor const &color) {
    Vec2f bboxmin( std::numeric_limits<float>::max(),  std::numeric_limits<float>::max());
    Vec2f bboxmax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    Vec2f clamp(image.get_width()-1, image.get_height()-1);
    for (int i=0; i<3; i++) {
        for (int j=0; j<2; j++) {
            bboxmin[j] = std::max(0.f,      std::min(bboxmin[j], pts[i][j]));
            bboxmax[j] = std::min(clamp[j], std::max(bboxmax[j], pts[i][j]));
        }
    }
    Vec3f P;
    for (P.x=bboxmin.x; P.x<=bboxmax.x; ++P.x) {
        for (P.y=bboxmin.y; P.y<=bboxmax.y; ++P.y) {
            Vec3f bc_screen  = barycentric(pts[0], pts[1], pts[2], P);
            if (bc_screen.x<0 || bc_screen.y<0 || bc_screen.z<0) continue;
            P.z = 0;
            for (int i=0; i<3; ++i) {
                P.z += pts[i][2]*bc_screen[i];
            }
            auto idx = static_cast<size_t>(P.x + P.y * width);
            if(zbuffer[idx] < P.z) {
                zbuffer[idx] = P.z;
                //image.set(P.x, P.y, texture.get(u*texture.get_width(), v*texture.get_height()));
                image.set(P.x, P.y, color);
            }
        }
    }
}