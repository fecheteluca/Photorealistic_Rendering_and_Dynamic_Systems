#include "renderer.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

// Renderer class definitions
Renderer::Renderer(
    Scene aux_scene, 
    Camera aux_camera,
    int aux_nr_rays
) {
    scene = aux_scene;
    camera = aux_camera;
    W = camera.get_width();
    H = camera.get_height();

    // Create an image buffer with all pixels initialized to black
    image.resize(W * H * 3, 0);
    std::fill(image.begin(), image.end(), 0);

    nr_rays = aux_nr_rays;
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

void Renderer::display_image() {
    stbi_write_png("../../images/image.png", W, H, 3, image.data(), 0);
}

#pragma omp declare reduction(+: Vector : \
    omp_out=omp_out+omp_in) \
    initializer(omp_priv=Vector(0.0,0.0,0.0))

void Renderer::render() {
    int intensity_depth = INTENSITY_DEPTH;
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Ray ray = camera.get_ray(i, j);
            Intersection intersection = scene.get_closest_hit(ray);
            if (intersection.flag) {
                Vector vec_albedo_average = Vector();

                #pragma omp parallel for reduction(+:vec_albedo_average)
                for (int k = 0; k < nr_rays; k++) { 
                    Vector vec_albedo = scene.get_intensity(ray, intensity_depth);
                    vec_albedo_average = vec_albedo_average + vec_albedo;
                }
                
                vec_albedo_average = (1.0 / (nr_rays * 1.0)) * vec_albedo_average;
                set_image_pixel_color(i, j, vec_albedo_average);
            }
            
        }
    }
    display_image();
}
