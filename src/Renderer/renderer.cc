#include "renderer.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

// Renderer class definitions
Renderer::Renderer(Scene aux_scene, Camera aux_camera) {
    scene = aux_scene;
    camera = aux_camera;
    W = camera.get_width();
    H = camera.get_height();

    // Create an image buffer with all pixels initialized to black
    image.resize(W * H * 3, 0);
    std::fill(image.begin(), image.end(), 0);
}

void Renderer::set_image_pixel_color(int i, int j, Vector vec_albedo) {
    image[(i * W + j) * 3 + 0] = vec_albedo.get_x() * 255; // Red
    image[(i * W + j) * 3 + 1] = vec_albedo.get_y() * 255; // Green
    image[(i * W + j) * 3 + 2] = vec_albedo.get_z() * 255; // Blue
}

Vector Renderer::get_image_pixel_color(int i, int j) {
    return Vector(image[(i * W + j) * 3 + 0], image[(i * W + j) * 3 + 1], image[(i * W + j) * 3 + 2]);
}

void Renderer::display_image() {
    stbi_write_png("../../images/image.png", W, H, 3, image.data(), 0);
}

void Renderer::render_default_image() {
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            set_image_pixel_color(i, j, Vector(0, 255, 0));
        }
    }
    display_image();
}

void Renderer::render_ray_scene_intersection() {
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Ray ray = camera.get_ray(i, j);
            Intersection intersection = scene.get_closest_hit(ray);
            set_image_pixel_color(i, j, intersection.vec_albedo);
        }
    }
    display_image();
}
