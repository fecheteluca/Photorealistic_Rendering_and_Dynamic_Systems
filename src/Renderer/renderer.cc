#include "renderer.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

// Renderer class definitions
Renderer::Renderer(Scene aux_scene, Camera aux_camera, Vector aux_vec_light_source, double aux_light_intensity) {
    scene = aux_scene;
    camera = aux_camera;
    W = camera.get_width();
    H = camera.get_height();
    vec_light_source = aux_vec_light_source;
    light_intensity = aux_light_intensity;

    // Create an image buffer with all pixels initialized to black
    image.resize(W * H * 3, 0);
    std::fill(image.begin(), image.end(), 0);
}

void Renderer::set_image_pixel_color(int i, int j, Vector vec_albedo) {
    double gamma = 2.2;
    image[(i * W + j) * 3 + 0] = pow(vec_albedo.get_x(), 1.0 / gamma) * 255; // Red
    image[(i * W + j) * 3 + 1] = pow(vec_albedo.get_y(), 1.0 / gamma) * 255; // Green
    image[(i * W + j) * 3 + 2] = pow(vec_albedo.get_z(), 1.0 / gamma) * 255; // Blue
}

Vector Renderer::get_image_pixel_color(int i, int j) {
    return Vector(image[(i * W + j) * 3 + 0], image[(i * W + j) * 3 + 1], image[(i * W + j) * 3 + 2]);
}

Vector Renderer::get_reflected_intensity(Intersection intersection) {
    double distance = (vec_light_source - intersection.vec_point).norm();
    Vector omega    = vec_light_source - intersection.vec_point;
    omega.normalize();

    double first_factor = light_intensity / (4 * M_PI * pow(distance, 2));
    Vector second_factor = (1 / M_PI) * intersection.vec_albedo;

    double third_factor = 0;
    Vector vec_offseted_point = intersection.vec_point + 0.001 * intersection.vec_normal;
    Ray ray = Ray(vec_offseted_point, omega);   
    Intersection src_pnt_intersection = scene.get_closest_hit(ray);
    if (distance <= src_pnt_intersection.distance) {
        third_factor = 1;
    }

    double fourth_factor = dot(intersection.vec_normal, omega);
    return first_factor * second_factor * third_factor * fourth_factor;
}

void Renderer::display_image() {
    stbi_write_png("../../images/image.png", W, H, 3, image.data(), 0);
}

// Render speific methods for displaying images
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
            if (intersection.flag) {
                set_image_pixel_color(i, j, intersection.vec_albedo);
            }
        }
    }
    display_image();
}

void Renderer::render_shading_and_shadows() {
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Ray ray = camera.get_ray(i, j);
            Intersection intersection = scene.get_closest_hit(ray);
            if (intersection.flag) {
                Vector vec_albedo_pixel = get_reflected_intensity(intersection);
                set_image_pixel_color(i, j, vec_albedo_pixel);
            }
        }
    }
    display_image();
}
