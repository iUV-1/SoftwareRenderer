//
// Created by iUV on 3/7/2025.
//
#include "my_gl.hpp"
#include "../../../Utilities/Adapters/TGAtoSDLAdapter.hpp"
#include "../../../Utilities/Interfaces/IShader.hpp"

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

// Overload with vec3
void line(Vec3f t0, Vec3f t1, TGAImage &image, TGAColor color) {
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

// Calculate barycentric value of a point in a triangle
// Returns (-1, 1, 1) in case the triangle is degenerate
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

    return Vec3f(-1,1,1); // in this case generate negative coordinates, it will be thrown away by the rasterizator
}

void triangle(Vec3f *pts, TGAImage &image, float *zbuffer, int width, IShader &shader) {
    Vec2f bboxmin( MAX_FLOAT,  MAX_FLOAT);
    Vec2f bboxmax(-MAX_FLOAT, -MAX_FLOAT);
    //std::cout << pts[0] << pts[1] << pts[2] << std::endl;
    Vec2f clamp(image.get_width()-1, image.get_height()-1);
    for (int i=0; i<3; i++) {
        for (int j=0; j<2; j++) {
            bboxmin[j] = std::max(0.f,      std::min(bboxmin[j], pts[i][j]));
            bboxmax[j] = std::min(clamp[j], std::max(bboxmax[j], pts[i][j]));
        }
    }
    //Vec3f P;
    for (int i = bboxmin.x; i <=bboxmax.x; i++) {
        for (int j = bboxmin.y; j<=bboxmax.y; j++) {
            Vec3f P(i, j, 0);
            Vec3f bc_screen  = barycentric(pts[0], pts[1], pts[2], P);
            if (bc_screen.x<0 || bc_screen.y<0 || bc_screen.z<0) continue;
            P.z = 0;
            for (int i=0; i<3; ++i) {
                P.z += pts[i][2]*bc_screen[i];
            }
            auto idx = static_cast<size_t>(P.x + P.y * width);
            if(zbuffer[idx] < P.z) {
                zbuffer[idx] = P.z;
                // Use shader
                TGAColor color;
                shader.fragment(bc_screen, color);
                image.set(P.x, P.y, color);
            }
        }
    }
}

void triangle(Vec3f *pts, SDL_Surface *surface, float *zbuffer, int width, IShader &shader) {
    Vec2f bboxmin( MAX_FLOAT,  MAX_FLOAT);
    Vec2f bboxmax(-MAX_FLOAT, -MAX_FLOAT);
    //std::cout << pts[0] << pts[1] << pts[2] << std::endl;
    Vec2f clamp(surface->w -1, surface->h -1);
    for (int i=0; i<3; i++) {
        for (int j=0; j<2; j++) {
            bboxmin[j] = std::max(0.f,      std::min(bboxmin[j], pts[i][j]));
            bboxmax[j] = std::min(clamp[j], std::max(bboxmax[j], pts[i][j]));
        }
    }
    //Vec3f P;
    for (int i = bboxmin.x; i <=bboxmax.x; i++) {
        for (int j = bboxmin.y; j<=bboxmax.y; j++) {
            Vec3f P(i, j, 0);
            Vec3f bc_screen  = barycentric(pts[0], pts[1], pts[2], P);
            if (bc_screen.x<0 || bc_screen.y<0 || bc_screen.z<0) continue;
            P.z = 0;
            for (int i=0; i<3; ++i) {
                P.z += pts[i][2]*bc_screen[i];
            }
            auto idx = static_cast<size_t>(P.x + P.y * width);
            if(zbuffer[idx] < P.z) {
                zbuffer[idx] = P.z;
                // Use shader
                TGAColor color;
                shader.fragment(bc_screen, color);
                TGAtoSDLAdapter::SetTGAPixel(surface, P.x, P.y, color);
                //image.set(P.x, P.y, color);
            }
        }
    }
}

void triangle(Vec3f *pts, SDL_Surface *surface, RectBuffer &zbuffer, IShader &shader) {
    Vec2f bboxmin( MAX_FLOAT,  MAX_FLOAT);
    Vec2f bboxmax(-MAX_FLOAT, -MAX_FLOAT);
    //std::cout << pts[0] << pts[1] << pts[2] << std::endl;
    Vec2f clamp(surface->w -1, surface->h -1);
    for (int i=0; i<3; i++) {
        for (int j=0; j<2; j++) {
            bboxmin[j] = std::max(0.f,      std::min(bboxmin[j], pts[i][j]));
            bboxmax[j] = std::min(clamp[j], std::max(bboxmax[j], pts[i][j]));
        }
    }
    //Vec3f P;
    for (int i = bboxmin.x; i <=bboxmax.x; i++) {
        for (int j = bboxmin.y; j<=bboxmax.y; j++) {
            Vec3f P(i, j, 0);
            Vec3f bc_screen  = barycentric(pts[0], pts[1], pts[2], P);
            if (bc_screen.x<0 || bc_screen.y<0 || bc_screen.z<0) continue;
            P.z = 0;
            for (int i=0; i<3; ++i) {
                P.z += pts[i][2]*bc_screen[i];
            }
            auto idx = static_cast<size_t>(P.x + P.y * zbuffer.width);
            if(zbuffer[idx] < P.z) {
                zbuffer[idx] = P.z;
                // Use shader
                TGAColor color;
                shader.fragment(bc_screen, color);
                TGAtoSDLAdapter::SetTGAPixel(surface, P.x, P.y, color);
                //image.set(P.x, P.y, color);
            }
        }
    }
}

void triangle(Vec3f *pts, Window *window, RectBuffer &zbuffer, IShader &shader) {
    Vec2f bboxmin( MAX_FLOAT,  MAX_FLOAT);
    Vec2f bboxmax(-MAX_FLOAT, -MAX_FLOAT);

    Vec2f clamp(window->mWidth -1, window->mHeight -1);

    for (int i=0; i<3; i++) {
        for (int j=0; j<2; j++) {
            bboxmin[j] = std::max(0.f,      std::min(bboxmin[j], pts[i][j]));
            bboxmax[j] = std::min(clamp[j], std::max(bboxmax[j], pts[i][j]));
        }
    }

//#pragma omp parallel for schedule(guided, 64)
        for(int j = static_cast<int>(bboxmin.y); j <= static_cast<int>(bboxmax.y); j++) {
            for(int i = static_cast<int>(bboxmin.x); i <= static_cast<int>(bboxmax.x); i++) {

            Vec3f P(i, j, 0);
            Vec3f bc_screen  = barycentric(pts[0], pts[1], pts[2], P);
            if (bc_screen.x<0 || bc_screen.y<0 || bc_screen.z<0) continue;
            P.z = 0;
            for (int k=0; k<3; ++k) {
                P.z += pts[k].z * bc_screen[k];
            }
            auto idx = static_cast<size_t>(P.x + P.y * zbuffer.width);
            if(zbuffer[idx] < P.z) {
//#pragma omp atomic write
                zbuffer[idx] = P.z;
                // Use shader
                TGAColor color;
                shader.fragment(bc_screen, color);
                window->Plot(P.x, P.y, color.r, color.g, color.b, color.a);
            }
        }
    }
}

// Only for depth buffer
void triangle(Vec3f *pts, int w, int h, RectBuffer &zbuffer) {
    Vec2f bboxmin( MAX_FLOAT,  MAX_FLOAT);
    Vec2f bboxmax(-MAX_FLOAT, -MAX_FLOAT);

    Vec2f clamp(w -1, h -1);

    for (int i=0; i<3; i++) {
        for (int j=0; j<2; j++) {
            bboxmin[j] = std::max(0.f,      std::min(bboxmin[j], pts[i][j]));
            bboxmax[j] = std::min(clamp[j], std::max(bboxmax[j], pts[i][j]));
        }
    }

//#pragma omp parallel for schedule(guided, 64)
    for(int i = static_cast<int>(bboxmin.x); i <= static_cast<int>(bboxmax.x); i++) {
        for(int j = static_cast<int>(bboxmin.y); j <= static_cast<int>(bboxmax.y); j++) {
            Vec3f P(i, j, 0);
            Vec3f bc_screen  = barycentric(pts[0], pts[1], pts[2], P);
            if (bc_screen.x<0 || bc_screen.y<0 || bc_screen.z<0) continue;
            P.z = 0;
            for (int i=0; i<3; ++i) {
                P.z += pts[i].z*bc_screen[i];
            }
            auto idx = static_cast<size_t>(P.x + P.y * zbuffer.width);
            if(zbuffer[idx] < P.z) {
//#pragma omp atomic write
                zbuffer[idx] = P.z;
            }
        }
    }
}

// Draw a triangle in wireframe mode
void wireframe_trig(Vec3f *pts, TGAImage &image, TGAColor color) {
    line(pts[0] ,pts[1], image, color);
    line(pts[1], pts[2], image, color);
    line(pts[0], pts[2], image, color);
}