//
// Created by iUV on 3/7/2025.
//

#include "my_gl.hpp"
#include "geometry.h"
#include "tgaimage.h"

Matrix4x4f ModelView;
Matrix4x4f Projection;
Matrix4x4f Viewport;

// Create a homogonized matrix from a vector
Matrix<float> homogonize(Vec3f v) {
    Matrix<float> result(4, 1);
    result[0][0] = v.x;
    result[1][0] = v.y;
    result[2][0] = v.z;
    result[3][0] = 1;
    return result;
}

IShader::~IShader() {

}


// De-homogonize it
Vec3f dehomogonize(Matrix<float> const &m) {
    Vec3f result;
    result.x = m[0][0] / m[3][0];
    result.y = m[1][0] / m[3][0];
    result.z = m[2][0] / m[3][0];
    return result;
}

// Old function. Meant to project the points using a projection matrix
Vec3f project(Vec3f v) {
    // row matrix
    Matrix<float> homogonized = homogonize(v);
    //Matrix<float> transformed = transfrom.multiply(homogonized); // Matrix4x4f multiply by homogonized vector
    Matrix<float> transformed = Projection*homogonized;
    Vec3f result = dehomogonize(transformed);
    return result;
}

void Project(float coeff) {
    Projection = Matrix4x4f::identity();
    Projection[3][2] = -1/coeff;
}

// Similar to gluLookAt, create a camera transformation matrix
// Formula (8.4) in textbook
void LookAt(Vec3f eye, Vec3f center, Vec3f up) {
    Vec3f z = (eye-center).normalize();
    Vec3f x = (up^z).normalize();
    Vec3f y = (z^x).normalize();
    // Inverse rotation matrix
    auto Minv = Matrix4x4f::identity();
    // Translation matrix
    auto Tr = Matrix4x4f::identity();
    for (int i = 0; i < 3; i++) {
        Minv[0][i] = x[i];
        Minv[1][i] = y[i];
        Minv[2][i] = z[i];
        Tr[i][3] = -eye[i];
    }
    ModelView = Minv * Tr;
    //return M_cam;
}

// Matrix representation of viewport transformation
// Also includes depth because viewport is a box
void world2screen(Vec3f v, int w, int h, float depth) {
    Viewport = Matrix4x4f::identity();
    Viewport[0][3] = v.x+w/2.f;
    Viewport[1][3] = v.y+h/2.f;
    Viewport[2][3] = depth/2.f;

    Viewport[0][0] = w/2.f;
    Viewport[1][1] = h/2.f;
    Viewport[2][2] = depth/2.f;
}

// Set Viewport Matrix (M_vp)
// Section 8.1 in textbook
void SetViewport(int width, int height, float depth) {
    Viewport = Matrix4x4f::identity();
    Viewport[0][0] = width/2.f;
    Viewport[1][1] = height/2.f;
    Viewport[0][3] = (width-1)/2.f;
    Viewport[1][3] = (height-1)/2.f;

    Viewport[2][3] = depth/2.f;
    Viewport[2][2] = depth/2.f;
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
                // Use shader
                TGAColor color;
                shader.fragment(bc_screen, color);
                image.set(P.x, P.y, color);
            }
        }
    }
}


// Overload that is just intensity
void triangle(Vec3f *pts, TGAImage &image, float *zbuffer, TGAColor const &color, int width) {
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

/* retirement home */

/*
Vec3f world2screen(Vec3f v) {
    return Vec3f(static_cast<int>((v.x+1.)*width/2.+.5), static_cast<int>((v.y+1.)*height/2. + .5), v.z);
}*/
