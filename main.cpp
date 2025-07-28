#include <vector>
#include <iostream>
#include <limits>
#include <chrono>
#include <ctime>
#include <sstream>

#include "model.h"
#include "geometry.h"
#include "my_gl.hpp"
#include "tgaimage.h"

float c;

Model *model = nullptr;
TGAImage tex_file(1024 ,1024,TGAImage::RGB);
TGAImage normal_file(1024, 1024, TGAImage::RGB);
TGAImage specular_file(1024, 1024, TGAImage::RGB);
bool use_normal_map = false;
bool use_specular_map = false;

const int width = 800;
const int height = 800;

const float depth = 100.0f; // Far clipping plane

Vec3f light = Vec3f(-1.0, 1.0, 1.0).normalize();

struct GouraudShaderReference: IShader {
    Vec3f varying_intensity; // intensity of a vertex
    Matrix<float> varying_uv = Matrix<float>(2, 3); // 2x3 matrix containing uv coordinate of 3 vertex (a trig)
    Matrix4x4f uniform_M; // Projection*ModelView
    Matrix4x4f uniform_MIT; // same as above but invert_transpose()

    Matrix<float> vertex(int iface, int nthvert) override{
        //light.normalize();
        Vec3f v = model->vert(iface, nthvert);
        Vec3f n = model->normal(iface, nthvert);
        // Set the column of varying_uv to texture position in Vec2f
        varying_uv.set_col(nthvert, model->texcoord(iface, nthvert));
        // Cap at 0
        varying_intensity[nthvert] = std::max(0.f, n*light);
        return Viewport*uniform_M*homogonize(v, 1.);
    }
    // bar is the barycentric of that vertex
    bool fragment(Vec3f bar, TGAColor &color) override {
        // Cap at 0
        float intensity = std::max(0.f, varying_intensity*bar);
        // Convert barycentric vector to a matrix
        // NOTE: Somehow making a new variable is faster than making it inline?
        Matrix<float> bary = Matrix(bar); // 1x3 row matrix that represent a vector
        // Matrix<float> uv = varying_uv*Matrix(bar) <- Slower!
        Matrix<float> uv = varying_uv*bary; // 1x2 Matrix (Basically a Vec2f)
        TGAColor texColor = tex_file.get(uv[0][0] * tex_file.get_width(), uv[1][0] * tex_file.get_height());
        color = texColor * intensity;
        return false;
    }
};

// Gouraud Shader uses vertex data to calculate light value
// However, this shader uses interpolation to compute the normal vector per pixel
// This is also called smooth shading
struct GouraudShader: IShader {
    Vec3f varying_intensity; // intensity of a vertex
    Matrix<float> varying_uv = Matrix<float>(2, 3); // 2x3 matrix containing uv coordinate of 3 vertex (a trig)
    Matrix4x4f uniform_M; // Projection*ModelView
    Matrix4x4f uniform_MIT; // same as above but invert_transpose()

    Matrix<float> vertex(int iface, int nthvert) override{
        //light.normalize();
        Vec3f v = model->vert(iface, nthvert);
        Vec3f n = model->normal(iface, nthvert);
        // Set the column of varying_uv to texture position in Vec2f
        varying_uv[0][nthvert] = model->texcoord(iface, nthvert).x;
        varying_uv[1][nthvert] = model->texcoord(iface, nthvert).y;
        // Cap at 0
        varying_intensity[nthvert] = std::max(0.f, n*light);
        return Viewport*Projection*ModelView*homogonize(v, 1.);
    }
    // bar is the barycentric of that vertex
    bool fragment(Vec3f bar, TGAColor &color) override {
        // Cap at 0
        // Convert barycentric vector to a matrix
        // NOTE: Somehow making a new variable is faster than making it inline?
        Matrix<float> bary = Matrix(bar); // 1x3 row matrix that represent a vector
        // Matrix<float> uv = varying_uv*Matrix(bar) <- Slower!
        Matrix<float> uv = varying_uv*bary; // 1x2 Matrix (Basically a Vec2f)

        Vec3f norm;
        if(false)
            norm = model->normal(uv[0][0], uv[1][0]);
        else {
            TGAColor normal_color = normal_file.get(uv[0][0] * normal_file.get_width(), uv[1][0] * normal_file.get_height());
            norm = Vec3f(normal_color.r, normal_color.g, normal_color.b);
        }
        // Transforming
        Vec3f l = dehomogonize(uniform_M  *homogonize(light, 0.f)).normalize();
        //l.z = -l.z;
        l.x = -l.x;
        l.y = -l.y;
        l.z = -l.z;
        Vec3f n = dehomogonize(uniform_MIT*homogonize(norm, 0.f)).normalize();
/*        if(varying_intensity*bar != n*l) {
            std::cout << "mismatched value" << std::endl;
        }*/
        float diff = std::max(0.f, n * l); // diffuse intensity value

        TGAColor texColor = tex_file.get(uv[0][0] * tex_file.get_width(), uv[1][0] * tex_file.get_height());
        color = texColor * diff;
        return false;
    }
};


// Phong Shader includes Phong reflection model and specular mapping
struct PhongShader: IShader {
    Matrix<float> varying_uv = Matrix<float>(2, 3); // 2x3 matrix containing uv coordinate of 3 vertex (a trig)
    Matrix4x4f uniform_M; // Projection*ModelView
    Matrix4x4f uniform_MIT; // same as above but invert_transpose()

    Matrix<float> vertex(int iface, int nthvert) override{
        Vec3f v = model->vert(iface, nthvert);
        // Set the column of varying_uv to texture position in Vec2f
        varying_uv.set_col(nthvert, model->texcoord(iface, nthvert));

        return Viewport*uniform_M*homogonize(v, 1.);
    }
    // bar is the barycentric of that vertex
    bool fragment(Vec3f bar, TGAColor &color) override {
        // Convert barycentric vector to a matrix
        // NOTE: Somehow making a new variable is faster than making it inline?
        Matrix<float> bary = Matrix(bar); // 1x3 row matrix that represent a vector
        // Matrix<float> uv = varying_uv*Matrix(bar) <- Slower!
        Matrix<float> mat_uv = varying_uv*bary; // 1x2 Matrix (Basically a Vec2f)
        Vec2f uv = Vec2f(mat_uv[0][0], mat_uv[1][0]);
        // Get the normal vector of that mesh based on the setting
        Vec3f norm;
        if(!use_normal_map)
            norm = model->normal(uv.u, uv.v);
        else {
            TGAColor normal_color = normal_file.get(uv.u * normal_file.get_width(), uv.v * normal_file.get_height());
            norm = Vec3f(normal_color.r, normal_color.g, normal_color.b);
        }
        /// Insanely costly calculations
        // Transform the normal vector to the eye space
        Vec3f n = dehomogonize(uniform_MIT*homogonize(norm, 0.f)).normalize();
        // Same as above
        Vec3f l = dehomogonize(uniform_M  *homogonize(light, 0.f)).normalize();
        l.z = -l.z; // I think this formula is meant to work for a different axis?
        l.y = -l.y;
        l.x = -l.x;
        Vec3f r = (n*(n*l*2.f) - l).normalize(); // reflection vector
        float diff = std::max(0.f, n * l); // diffuse intensity value
        // Specular
        float spec = 0.f;
        if(use_specular_map) {
            float spec_map_val = specular_file.get(uv.u * specular_file.get_width(), uv.v * specular_file.get_height()).r;
            spec = pow(std::max(r.z, 0.f), spec_map_val);
        } else {
            spec = std::max(r.z, 0.f);
        }
        TGAColor texColor = tex_file.get(uv.u * tex_file.get_width(), uv.v * tex_file.get_height());
        for (int i = 0; i < 3; i++) color[i] = std::min<float>(5 + texColor[i]*(diff + .6*spec), 255.f);
        color = texColor * diff;

        return false;
    }
};


// Like above but includes a shadow pass
struct PhongShaderShadow: IShader {
    Matrix<float> varying_uv = Matrix<float>(2, 3); // 2x3 matrix containing uv coordinate of 3 vertex (a trig)
    Matrix3x3f varying_tri; // 3x3 matrix containing verticies of a trig
    float* depth_buffer;
    Matrix4x4f uniform_M; // Projection*ModelView
    Matrix4x4f uniform_MIT; // same as above but invert_transpose()
    Matrix4x4f uniform_Mshadow; // Shadow transformation

    Matrix<float> vertex(int iface, int nthvert) override{
        Vec3f v = model->vert(iface, nthvert);
        // Set the column of varying_uv to texture position in Vec2f
        varying_uv.set_col(nthvert, model->texcoord(iface, nthvert));
        // Set the column of vertex the triangle using vert index
        Matrix<float> transformed_vert = Viewport*uniform_M*homogonize(v, 1.);
        varying_tri.set_col(nthvert, dehomogonize(transformed_vert));

        return transformed_vert;
    }
    // bar is the barycentric of that vertex
    bool fragment(Vec3f bar, TGAColor &color) override {
        // Convert barycentric vector to a matrix
        // NOTE: Somehow making a new variable is faster than making it inline?
        Matrix<float> bary = Matrix(bar); // 1x3 row matrix that represent a vector
        // Matrix<float> uv = varying_uv*Matrix(bar) <- Slower!
        Matrix<float> mat_uv = varying_uv*bary; // 1x2 Matrix (Basically a Vec2f)

        // Get shadow position from buffer
        Matrix<float> matrix_p = varying_tri * bary;
        Vec3f p = { matrix_p[0][0], matrix_p[1][0], matrix_p[2][0]};
        Matrix<float> shadow_buffer_pt = uniform_Mshadow* homogonize(p, 1.f);
        Vec3f shadow_p = dehomogonize(shadow_buffer_pt);
        auto shadow_buf_idx = static_cast<size_t>(shadow_p.x + shadow_p.y * width);
        float shadow = .3 + 7*(depth_buffer[shadow_buf_idx] < shadow_p.z);

        Vec2f uv = Vec2f(mat_uv[0][0], mat_uv[1][0]);
        // Get the normal vector of that mesh based on the setting
        Vec3f norm;
        if(!use_normal_map)
            norm = model->normal(uv.u, uv.v);
        else {
            TGAColor normal_color = normal_file.get(uv.u * normal_file.get_width(), uv.v * normal_file.get_height());
            norm = Vec3f(normal_color.r, normal_color.g, normal_color.b);
        }
        /// Insanely costly calculations
        // Transform the normal vector to the eye space
        Vec3f n = dehomogonize(uniform_MIT*homogonize(norm, 0.f)).normalize();
        // Same as above
        Vec3f l = dehomogonize(uniform_M  *homogonize(light, 0.f)).normalize();
        l.z = -l.z; // I think this formula is meant to work for a different axis?
        l.y = -l.y;
        l.x = -l.x;
        Vec3f r = (n*(n*l*2.f) - l).normalize(); // reflection vector
        float diff = std::max(0.f, n * l); // diffuse intensity value
        // Specular
        float spec = 0.f;
        if(use_specular_map) {
            float spec_map_val = specular_file.get(uv.u * specular_file.get_width(), uv.v * specular_file.get_height()).r;
            spec = pow(std::max(r.z, 0.f), spec_map_val);
        } else {
            spec = std::max(r.z, 0.f);
        }
        TGAColor texColor = tex_file.get(uv.u * tex_file.get_width(), uv.v * tex_file.get_height());
        for (int i = 0; i < 3; i++) color[i] = std::min<float>(5 + texColor[i]*shadow*(1.2*diff + .6*spec), 255.f);
        color = texColor * diff;

        return false;
    }
};

// The famous Rainbow Triangle
// vertex shader is discarded entirely
struct RainbowShader: IShader {
    Matrix<float> vertex(int iface, int nthvert) override{
        Matrix<float> dummy = Matrix4x4f();
        return dummy;
    }

    bool fragment(Vec3f bar, TGAColor &color) override{
        TGAColor rainbow(bar.x * 255, bar.y * 255, bar.z * 255, 255);
        color = rainbow;
        return false;
    }
};

// Copy zbuffer to a framebuffer (Image in this case)
struct DepthShader: IShader {
    Matrix4x4f uniform_M; // Projection*ModelView
    Matrix3x3f varying_tri; // 3x3 matrix containing vertex position of a trig

    // Typical vertex rendering
    Matrix<float> vertex(int iface, int nthvert) override{
        Vec3f v = model->vert(iface, nthvert);
        // Set the column of varying_uv to texture position in Vec2f
        Matrix<float> transformed_vert = Viewport*uniform_M*homogonize(v, 1.);
        varying_tri.set_col(nthvert, dehomogonize(transformed_vert));
        return transformed_vert;
    }

    //
    bool fragment(Vec3f bar, TGAColor &color) override {
        Matrix<float> bary(bar);
        Matrix<float> pt = varying_tri * bary; // point interpolated after transformation

        // Set the brightness based on how far is it from the camera
        // clamp
        float dist = std::clamp(pt[2][0]/depth, 0.f, 1.f);
        color = TGAColor(255, 255, 255) *
                (dist);

        return false;
    }
};

Vec3f rasterize(IShader *shader, int iface, int nthvert) {
    // Apply vertex shader
    Matrix<float> homogonized = shader->vertex(iface, nthvert);
    Vec3f result = dehomogonize(homogonized);
    // Round the result to apply to screen
    result.x = std::round(result.x);
    result.y = std::round(result.y);
    result.z = std::round(result.z);

    return result;
}

int main(int argc, char** argv) {

    /* the famous rainbow triangle */
    if (false) {
        float *zbuffer = new float[width*height];
        std::fill(zbuffer, zbuffer + width*height, -std::numeric_limits<float>::max()); // set every value in zbuffer to -inf

        Vec3f pt1(0, 0, 0);
        Vec3f pt2(128, 256, 0);
        Vec3f pt3(256, 0, 0);
        Vec3f trig[3] = {pt1, pt2, pt3};

        RainbowShader shader = RainbowShader();

        auto frame = new TGAImage(width, width, TGAImage::RGB);

        triangle(trig, *frame, zbuffer, width, shader);
        frame->flip_vertically();
        // Get timing of the render
        frame->write_tga_file("rainbow triangle.tga");
        delete[] zbuffer;
        delete frame;
        return 0;
    }

    // Time the render
    auto before = std::chrono::system_clock::now();

    if(argc < 2)
        model = new Model("obj/african_head.obj");
    else
        model = new Model(argv[1]);


    //TGAImage tex_file(1024,1024,TGAImage::RGB);
    if(argc < 3)
        tex_file.read_tga_file("obj/UV Grid.tga");
    else
        tex_file.read_tga_file(argv[2]);


    tex_file.flip_vertically();

    // Check if normal map is included in the args
    // If not, use the model embeded normal
    if(argc >= 3) {
        normal_file.read_tga_file(argv[3]);
        normal_file.flip_vertically();
        use_normal_map = true;
    }

    if(argc >= 4) {
        specular_file.read_tga_file(argv[4]);
        specular_file.flip_vertically();
        use_specular_map = true;
    }

    // camera setting
    Vec3f eye(0, 0, 2);
    Vec3f cam(0, 0, 0);
    Vec3f up(0, 1, 0);

    /*Vec3f eye(3, 2, 3);
    Vec3f cam(0, 0, 0);
    Vec3f up(0, 1, 0);*/

    auto render = std::chrono::system_clock::now();

    /* Depth map*/
    // Init depth buffer
    auto *depth_buffer_arr = new float[width*height];
    std::fill(depth_buffer_arr, depth_buffer_arr + width*height, -std::numeric_limits<float>::max()); // set every value in zbuffer to -inf

    auto depth_buffer = TGAImage(width, height, TGAImage::RGB);
    // Init shader
    LookAt(light, cam, up); // Render from the light (normalized)
    Project(0); // Render light in orthographic mode
    SetViewport(width, height, 255.0f);

    SetViewport(width / 8, height/8, width * 3/4, height * 3/4, depth); // Clamp the image into the center with margins (3/4 of the screen)

    DepthShader depth_shader = DepthShader();
    depth_shader.uniform_M = Projection*ModelView;
    // Render
    for(int i = 0; i < model->nfaces(); ++i) {
        Vec3f screen_coords[3];
        for (int j=0; j<3; ++j)
            screen_coords[j] = rasterize(&depth_shader, i, j);

        triangle(screen_coords, depth_buffer, depth_buffer_arr, width, depth_shader);
    }
    depth_buffer.flip_vertically();
    Matrix4x4f M_Shadow = Viewport*depth_shader.uniform_M;

    /* Render */
    auto frame = TGAImage(width, height, TGAImage::RGB);

    // Setup GL
    LookAt(eye, cam, up);
    Project(-1/(eye-cam).norm());
    SetViewport(width, height, 255.0f);

    // Setup zbuffer
    float *zbuffer = new float[width*height];
    std::fill(zbuffer, zbuffer + width*height, -std::numeric_limits<float>::max()); // set every value in zbuffer to -inf

    // GouraudShader shader = GouraudShader();
    // PhongShader shader = PhongShader();
    // GouraudShaderReference shader = GouraudShaderReference();
    PhongShaderShadow shader = PhongShaderShadow();
    shader.uniform_M = Projection*ModelView;
    shader.uniform_MIT = shader.uniform_M;
    shader.uniform_MIT.inverseTranspose();
    shader.uniform_Mshadow = M_Shadow;
    shader.depth_buffer = depth_buffer_arr;

    for (int i=0; i<model->nfaces(); ++i) {
        Vec3f screen_coords[3];

        for (int j=0; j<3; ++j)
            screen_coords[j] = rasterize(&shader, i, j);

        // calculate normal
        // ^ is an overloaded operator that performs cross product calculation
        // world_coords[2] - world_coords[0] and the other are 2 vectors pointing from point
        // world_coords[0].
        Vec3f n = (screen_coords[2]-screen_coords[0])^(screen_coords[1]-screen_coords[0]);
        n.normalize();
        // calculate eye intensity by dot product between normal and eye vector
        float view_dir_intensity = eye*n;
        // back face culling

        if (view_dir_intensity<1) {
            triangle(screen_coords, frame, zbuffer, width, shader);
            //wireframe_trig(screen_coords, frame, TGAColor(255, 255, 255, 255));
        }
    }
    // set origin to the bottom left corner
    frame.flip_vertically();

    // Get timing of the render
    auto now = std::chrono::system_clock::now();
    auto finish_time = std::chrono::system_clock::to_time_t(now);
    std::tm local_time = *std::localtime(&finish_time);
    std::stringstream sstream;
    sstream << "../output/" << local_time.tm_mon + 1 << "-" << local_time.tm_mday << "_"
            << local_time.tm_hour << "-" << local_time.tm_min << ".tga";
    frame.write_tga_file(sstream.str().c_str());
    sstream << "_buffer.tga";
    depth_buffer.write_tga_file(sstream.str().c_str());

    // How long the render takes
    std::chrono::duration<double> diff = now - before;
    auto ms = duration_cast<std::chrono::milliseconds>(diff);
    std::cout << "Elapsed time: " << ms.count() << " ms" << std::endl;
    diff = now-render;
    ms = duration_cast<std::chrono::milliseconds>(diff);
    std::cout << "Render time: " << ms.count() << " ms" << std::endl;

    delete model;
    delete[] zbuffer;
    delete[] depth_buffer_arr;
    return 0;
}