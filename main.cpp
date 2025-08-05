#include <iostream>
#include <chrono>
#include <ctime>
#include <sstream>

#include "model.h"
#include "geometry.h"
#include "my_gl.hpp"
#include "tgaimage.h"
#include "shaders.hpp"

constexpr int width = 800;
constexpr int height = 800;
constexpr float depth = 255.f;

Model *model;
TGAImage normal_file;
TGAImage tex_file;
TGAImage specular_file;

bool use_specular;
bool use_normal;

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
#if false
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
#endif

    // Time the render
    auto before = std::chrono::system_clock::now();

    if(argc < 2)
        model = new Model("obj/african_head.obj");
    else
        model = new Model(argv[1]);

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
        use_normal = true;
    }

    if(argc >= 4) {
        specular_file.read_tga_file(argv[4]);
        specular_file.flip_vertically();
        use_specular = true;
    }

    // camera setting
    Vec3f eye(0, 0, 2);
    Vec3f cam(0, 0, 0);
    Vec3f up(0, 1, 0);

    /*Vec3f eye(1, 1, 3);
    Vec3f cam(0, 0, 0);
    Vec3f up(0, 1, 0);*/

    /*Vec3f eye(3, 2, 3);
    Vec3f cam(0, 0, 0);
    Vec3f up(0, 1, 0);*/

    auto render = std::chrono::system_clock::now();

    /* Depth map */
    // Init depth buffer
    float *depth_buffer_arr = create_buffer(width, height);
    auto depth_buffer = TGAImage(width, height, TGAImage::RGB);
    // Init shader
    LookAt(light, cam, up); // Render from the light
    Project(0); // Render light in orthographic mode
    SetViewport(width / 8, height/8, width * 3./4, height * 3./4, depth); // Clamp the image into the center with margins (3/4 of the screen)

    DepthShader depth_shader = DepthShader();
    // Render
    for(int i = 0; i < model->nfaces(); ++i) {
        Vec3f screen_coords[3];
        for (int j=0; j<3; ++j)
            screen_coords[j] = rasterize(&depth_shader, i, j);

        Vec3f n = (screen_coords[2]-screen_coords[0])^(screen_coords[1]-screen_coords[0]);
        n.normalize();
        float view_dir_intensity = eye*n;

        if (view_dir_intensity<1)
            triangle(screen_coords, depth_buffer, depth_buffer_arr, width, depth_shader);
    }
    //depth_buffer.flip_vertically();
    Matrix4x4f M_Shadow = Viewport*Projection*ModelView;

    /* SSAO */
    /// Get depth from Camera
    auto depth_from_cam = TGAImage(width, height, TGAImage::RGB);
    // Setup GL
    LookAt(eye, cam, up);
    Project(-1/(eye-cam).norm());
    SetViewport(width/8, height/8, width*3./4, height*3./4, depth);

    // Setup zbuffer
    float *depth_cam_buf = create_buffer(width, height);
    depth_shader = DepthShader();
    for(int i = 0; i < model->nfaces(); ++i) {
        Vec3f screen_coords[3];
        for (int j=0; j<3; ++j)
            screen_coords[j] = rasterize(&depth_shader, i, j);

        Vec3f n = (screen_coords[2]-screen_coords[0])^(screen_coords[1]-screen_coords[0]);
        n.normalize();
        float view_dir_intensity = eye*n;

        if (view_dir_intensity<1)
            triangle(screen_coords, depth_from_cam, depth_cam_buf, width, depth_shader);
    }
    depth_from_cam.flip_vertically();
    /// SSAO
    auto ssao_pass = TGAImage(width, height, TGAImage::RGB);
    // Setup zbuffer
    auto ssao_buf = create_buffer(width, height);
    SSAOShader ssao_shader = SSAOShader();
    ssao_shader.depth = depth_from_cam;
    for(int i = 0; i < model->nfaces(); ++i) {
        Vec3f screen_coords[3];
        for (int j=0; j<3; ++j)
            screen_coords[j] = rasterize(&ssao_shader, i, j);

        Vec3f n = (screen_coords[2]-screen_coords[0])^(screen_coords[1]-screen_coords[0]);
        n.normalize();
        float view_dir_intensity = eye*n;

        if (view_dir_intensity<1)
            triangle(screen_coords, ssao_pass, ssao_buf, width, ssao_shader);
    }
    ssao_pass.flip_vertically();

    /* Render */
    auto frame = TGAImage(width, height, TGAImage::RGB);

    // Setup GL
    LookAt(eye, cam, up);
    Project(-1/(eye-cam).norm());
    SetViewport(width/8, height/8, width*3./4, height*3./4, depth);

    // Setup zbuffer
    float *zbuffer = create_buffer(width, height);
    // GouraudShader shader = GouraudShader();
    // PhongShader shader = PhongShader();
    // GouraudShaderReference shader = GouraudShaderReference();
    Matrix4x4f MVP = Viewport*Projection*ModelView;
    MVP.invert();
    // PhongShaderShadow shader = PhongShaderShadow(M_Shadow*MVP, depth_buffer_arr);
    PhongShader shader = PhongShader();

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
    sstream << "_ssao.tga";
    ssao_pass.write_tga_file(sstream.str().c_str());

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