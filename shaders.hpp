//
// Created by iuv on 7/30/25.
//

#ifndef SHADERS_HPP
#define SHADERS_HPP
//
// Created by iuv on 7/30/25.
//

#include <iostream>
#include <random>

#include "model.h"
#include "shaders.hpp"
#include "geometry.h"
#include "my_gl.hpp"

Vec3f light = Vec3f(-1.0, 1.0, 1.0).normalize();
extern Model *model;

extern TGAImage tex_file;
extern TGAImage normal_file;
extern TGAImage specular_file;

extern bool use_specular; // Use specular map
extern bool use_normal; // Use normal map

extern const float depth;
extern const int width;
extern const int height;

struct GouraudShaderReference: IShader {
    Vec3f varying_intensity; // intensity of a vertex

    Matrix<float> vertex(int iface, int nthvert) override{
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
        if(!use_normal)
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
        if(use_specular) {
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
    Matrix3x3<float> varying_tri; // 3x3 matrix containing verticies of a trig
    Matrix4x4f uniform_Mshadow; // Shadow transformation
    float* depth_buffer;

    PhongShaderShadow(Matrix4x4f uniform_shadow, float* depth_buffer) : uniform_Mshadow(std::move(uniform_shadow)), depth_buffer(depth_buffer) {}

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
        Vec2f uv = Vec2f(mat_uv[0][0], mat_uv[1][0]);

        // Get shadow position from buffer
        Matrix<float> matrix_p = varying_tri * bary;
        Vec3f p = { matrix_p[0][0], matrix_p[1][0], matrix_p[2][0]};
        Matrix<float> shadow_buffer_pt = uniform_Mshadow* homogonize(p, 1.f);
        Vec3f shadow_p = dehomogonize(shadow_buffer_pt);
        auto shadow_buf_idx = static_cast<size_t>(shadow_p.x + shadow_p.y * width);
        float shadow = .3 + 7*(depth_buffer[shadow_buf_idx] < shadow_p.z+43.34); // magic coeff to avoid z-fighting

        // Get the normal vector of that mesh based on the setting
        Vec3f norm;
        if(!use_normal)
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
        if(use_specular) {
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
    Matrix<float> vertex(int iface, int nthvert) override {
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
struct DepthShaderImage: IShader {
    Matrix3x3<float> varying_tri; // 3x3 matrix containing vertex position of a trig

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
        Matrix bary(bar);
        Matrix<float> pt = varying_tri * bary; // point interpolated after transformation

        // Set the brightness based on how far is it from the camera
        // clamp
        float dist = std::clamp(pt[2][0]/depth, 0.f, 1.f);
        color = TGAColor(255, 255, 255) *
                (dist);

        return false;
    }
};

struct DepthShader: IShader {
    // Typical vertex rendering
    Matrix<float> vertex(int iface, int nthvert) override{
        Vec3f v = model->vert(iface, nthvert);
        // Set the column of varying_uv to texture position in Vec2f
        // Matrix<float> transformed_vert = Viewport*uniform_M*homogonize(v, 1.);
        return Viewport*uniform_M*homogonize(v, 1.);
    }

    // Discard the fragment because we are only interested in the zbuffer produced by triangle()
    bool fragment(Vec3f bar, TGAColor &color) override {
        color = TGAColor(0,0,0,0);
        return false;
    }
};

constexpr int kernelSize = 16;
constexpr int noiseSize = 16;
constexpr float radius = 0.5f;
std::random_device rd;
std::mt19937 gen(rd());
/*
struct SSAOShader: IShader {
    Matrix3x3<float> varying_tri;
    Vec3f kernel[kernelSize];
    Vec3f noise[noiseSize];
    TGAImage depth;
    Vec2f uNoiseScale = Vec2f(width / 4, height / 4);
    TGAImage noise_texture = TGAImage(4,4,TGAImage::RGB);
    float *depth_buffer;

    SSAOShader() {
        // Generate random coordinates in the hemisphere
        for (int i = 0; i < kernelSize; ++i) {
            kernel[i] = Vec3f(
                random(-1.0f, 1.0f),
                random(-1.0f, 1.0f),
                random(0.0f, 1.0f)
                ).normalize();
            float scale = static_cast<float>(i) / static_cast<float>(kernelSize) * random(0.0f, 1.0f);
            scale = lerp(0.1f, 1.0f, scale * scale);
            kernel[i] *= scale;
        }
        // Generate rotational noise

        for(int i = 0; i < 4; ++i) {
            for(int j = 0; j < 4; ++j) {
                noise[i*4+j] = Vec3f(
                   random(-1.0f, 1.0f),
                   random(-1.0f, 1.0f),
                   0.0f
                ).normalize();
                noise_texture.set(i,j,TGAColor(
                        (noise[i*4+j].x * .5f + .5f) * 255,
                        (noise[i*4+j].y * .5f + .5f) * 255,
                        255));
            }
        }
    }

    Matrix3x3f varying_norm;
    // Typical vertex rendering
    Matrix<float> vertex(int iface, int nthvert) override{
        Vec3f v = model->vert(iface, nthvert);
        Vec3f n = model->normal(iface, nthvert);
        // Set the column of varying_uv to texture position in Vec2f
        varying_tri.set_col(nthvert, v);
        varying_norm.set_col(nthvert, v);
        return Viewport*uniform_M*homogonize(v, 1.);
    }

    float random(float r1, float r2) {
        std::uniform_real_distribution<float> dist(r1, r2);
        return dist(gen);
    }

    float lerp(float v0, float v1, float t) {
        return v0 + t * (v1 - v0);
    }


    bool fragment(Vec3f bar, TGAColor &color) override {
        // Get uv
        Matrix<float> bary = Matrix(bar); // 1x3 row matrix that represent a vector
        Matrix<float> mat_uv = varying_uv*bary; // 1x2 Matrix
        Vec2f uv = Vec2f(mat_uv[0][0], mat_uv[1][0]);

        // Get normal
        Vec3f norm;
        if(!use_normal)
            norm = model->normal(uv.u, uv.v);
        else {
            TGAColor normal_color = normal_file.get(uv.u * normal_file.get_width(), uv.v * normal_file.get_height());
            norm = Vec3f(normal_color.r, normal_color.g, normal_color.b);
        }

        // Get point
        Matrix<float> matrix_p = varying_tri * bary;
        Vec3f p = { matrix_p[0][0], matrix_p[1][0], matrix_p[2][0]};

        TGAColor rvac_color = noise_texture.get(uv.x * uNoiseScale.x,uv.y * uNoiseScale.y);
        Vec3f rvac = {
                static_cast<float>(rvac_color.r * 2.0 - 1.0),
                static_cast<float>(rvac_color.g * 2.0 - 1.0),
                0};
        Vec3f tangent = (rvac - norm * (rvac * norm));
        Vec3f bitangent = norm^tangent;
        // Create transformation matrix
        Matrix<float> M = Matrix<float>(3, 3);
        M.set_col(0, rvac);
        M.set_col(1, tangent);
        M.set_col(2, bitangent);

        /*
        // Transform kernel and test with depth buffer
        float occlusion = 0.f;
        for(int i = 0; i < kernelSize; i++) {
            // get sample position:
            Matrix<float> M_result = M * Matrix(kernel[i]);
            Vec3f sample = {M_result[0][0], M_result[1][0], M_result[2][0]};
            sample = sample * radius + p;

            // Project sample position
            Vec3f offset = dehomogonize(Projection*homogonize(sample, 1.f));
            offset.x = offset.x * 0.5f + 0.5f;
            offset.y = offset.y * 0.5f + 0.5f;

            // Get sample depth
            // Test with depth buffer
            //auto idx = static_cast<size_t>(offset.x + offset.y * width);
            //float sampleDepth = depth_buffer[idx];
            float sampleDepth = depth.get(offset.x * depth.get_width(), offset.y * depth.get_height()).r / 255.0f;
            if(sampleDepth > 0.) {
                std::cout << "something" << std::endl;
            }

            // range check & accumulate:
            float pDepth = p.z / 255.f;
            float rangeCheck= abs(pDepth - sampleDepth) < radius ? 1.0 : 0.0;
            if((sampleDepth <= sample.z ? 1.0 : 0.0) * rangeCheck >= 1) {
                //std::cout << "sampleDepth: " << sampleDepth << std::endl;
            }
            occlusion += (sampleDepth <= sample.z ? 1.0 : 0.0) * rangeCheck;
        }
        occlusion = 1.0 - (occlusion / kernelSize);

        // Produce image
        color = TGAColor(occlusion * 255, occlusion * 255, occlusion * 255);
        return false;

        float occlusion = 0.f;
        for (int i = 0; i < kernelSize; i++) {
            // get sample position
            Matrix<float> M_result = M * Matrix(kernel[i]);
            Vec3f sample = {M_result[0][0], M_result[1][0], M_result[2][0]};
            sample = sample * radius + p;

            // Project sample position to [0,1] screen coords
            Vec3f offset = dehomogonize(Projection * homogonize(sample, 1.f));
            offset.x = offset.x * 0.5f + 0.5f;
            offset.y = offset.y * 0.5f + 0.5f;

            // Fetch depth at this projected point
            size_t idx = static_cast<size_t>(offset.x * width + offset.y );
            if (idx > width*height) idx = width*height;
            float sampleDepth = depth_buffer[idx];

            // Normalize current sample.z the same way as depth shader
            //float sampleZNorm = ; // depth = far plane used in depth shader
            // Range check & accumulate
            float rangeCheck = fabs(p.z  - sampleDepth) < radius ? 1.0f : 0.0f;
            if (sampleDepth < sample.z) {
                occlusion += rangeCheck;
            }
        }
        occlusion = 1.0f - (occlusion / kernelSize);
        color = TGAColor(occlusion * 255, occlusion * 255, occlusion * 255);
        return false;
    }
};
*/


struct SSAOShader: IShader {
    Matrix<float> varying_pos_view = Matrix<float>(3, 3);
    Matrix<float> varying_norm_view = Matrix<float>(3, 3);
    Vec3f kernel[kernelSize];
    TGAImage depth;
    Vec2f uNoiseScale = Vec2f(width / 4, height / 4);
    TGAImage noise_texture = TGAImage(4,4,TGAImage::RGB);
    float *depth_buffer;

    float random(float r1, float r2) {
        std::uniform_real_distribution<float> dist(r1, r2);
        return dist(gen);
    }

    float lerp(float v0, float v1, float t) {
        return v0 + t * (v1 - v0);
    }
    SSAOShader() {
        // Generate random coordinates in the hemisphere
        for (int i = 0; i < kernelSize; ++i) {
            Vec3f sample(
                    random(-1.0f, 1.0f),
                    random(-1.0f, 1.0f),
                    random(0.0f, 1.0f)
            );
            sample.normalize();
            float scale = static_cast<float>(i) / static_cast<float>(kernelSize) * random(0.0f, 1.0f);
            scale = lerp(0.1f, 1.0f, scale * scale);
            kernel[i] = sample * scale;
        }
        // Generate rotational noise
        for(int i = 0; i < 4; ++i) {
            for(int j = 0; j < 4; ++j) {
                Vec3f n(
                        random(-1.0f, 1.0f),
                        random(-1.0f, 1.0f),
                        0.0f
                );
                n.normalize();
                noise_texture.set(i, j, TGAColor(
                        (n.x * .5f + .5f) * 255,
                        (n.y * .5f + .5f) * 255,
                        255)
                );
            }
        }
    }
    Matrix<float> vertex(int iface, int nthvert) override {
        Vec3f v = model->vert(iface, nthvert);
        Vec3f n = model->normal(iface, nthvert);

        // Store view-space position and normal for interpolation
        varying_pos_view.set_col(nthvert, dehomogonize(ModelView * homogonize(v, 1.f)));
        varying_norm_view.set_col(nthvert, dehomogonize(ModelView * homogonize(n, 0.f)).normalize());

        return Viewport * Projection * ModelView * homogonize(v, 1.f);
    }

    bool fragment(Vec3f bar, TGAColor &color) override {
        Matrix<float> bary = Matrix(bar);

        // Get interpolated view-space position and normal
        Matrix<float> frag_pos_matrix = varying_pos_view * bary;
        Vec3f frag_pos_view = {frag_pos_matrix[0][0], frag_pos_matrix[1][0], frag_pos_matrix[2][0]};

        Matrix<float> frag_norm_matrix = varying_norm_view * bary;
        Vec3f frag_norm_view = {frag_norm_matrix[0][0], frag_norm_matrix[1][0], frag_norm_matrix[2][0]};

        // Get noise vector
        TGAColor rvac_color = noise_texture.get(frag_pos_view.x * uNoiseScale.x, frag_pos_view.y * uNoiseScale.y);
        Vec3f rvac = {
                (rvac_color.r / 255.f) * 2.f - 1.f,
                (rvac_color.g / 255.f) * 2.f - 1.f,
                0.f
        };
        rvac.normalize();

        // Create TBN matrix for hemisphere rotation
        Vec3f tangent = (rvac - frag_norm_view * (rvac * frag_norm_view)).normalize();
        Vec3f bitangent = frag_norm_view ^ tangent;

        Matrix<float> TBN = Matrix<float>(3, 3);
        TBN.set_col(0, tangent);
        TBN.set_col(1, bitangent);
        TBN.set_col(2, frag_norm_view);

        float occlusion = 0.f;
        for (int i = 0; i < kernelSize; i++) {
            // Get sample position in view space
            Matrix<float> M_kernel = Matrix(kernel[i]);
            Matrix<float> M_sample_vec = TBN * M_kernel;
            Vec3f sample_vec = {M_sample_vec[0][0], M_sample_vec[1][0], M_sample_vec[2][0]};
            Vec3f sample_pos = frag_pos_view + sample_vec;

            // Project sample position to screen coords [0,1]
            Vec3f offset = dehomogonize(Projection * homogonize(sample_pos, 1.f));
            offset.x = offset.x * 0.5f + 0.5f;
            offset.y = offset.y * 0.5f + 0.5f;

            // Check if offset is within screen bounds
            if (offset.x < 0 || offset.y < 0 || offset.x > 1 || offset.y > 1) {
                continue;
            }

            // Fetch depth from depth buffer
            size_t idx = static_cast<size_t>(offset.x * width + offset.y);
            float sampleDepth = depth_buffer[idx];

            // Range check & accumulate
            float rangeCheck = fabs(frag_pos_view.z - sampleDepth) < radius ? 1.0f : 0.0f;
            if (sampleDepth > sample_pos.z) {
                occlusion += rangeCheck;
            }
        }
        occlusion = 1.0f - (occlusion / kernelSize);
        color = TGAColor(occlusion * 255, occlusion * 255, occlusion * 255);
        return false;
    }
};
#endif //SHADERS_HPP
