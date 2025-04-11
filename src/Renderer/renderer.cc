#include "renderer.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

// Renderer class definitions
Renderer::Renderer(
    Scene aux_scene, 
    Camera aux_camera,
    int aux_nr_rays,
    int aux_intensity_depth
) {
    scene = aux_scene;
    camera = aux_camera;
    W = camera.get_width();
    H = camera.get_height();

    // Create an image buffer with all pixels initialized to black
    image.resize(W * H * 3, 0);
    std::fill(image.begin(), image.end(), 0);

    nr_rays = aux_nr_rays;
    intensity_depth = aux_intensity_depth;
}

void Renderer::set_image_pixel_color(const int& i, const int& j, Vector& vec_albedo) {
    double gamma = 2.2;
    image[(i * W + j) * 3 + 0] = std::min(pow(vec_albedo.get_x(), 1.0 / gamma) * 254, 254.0); // Red
    image[(i * W + j) * 3 + 1] = std::min(pow(vec_albedo.get_y(), 1.0 / gamma) * 254, 254.0); // Green
    image[(i * W + j) * 3 + 2] = std::min(pow(vec_albedo.get_z(), 1.0 / gamma) * 254, 254.0); // Blue
}

Vector Renderer::get_image_pixel_color(const int& i, const int& j) {
    return Vector(image[(i * W + j) * 3 + 0], image[(i * W + j) * 3 + 1], image[(i * W + j) * 3 + 2]);
}

void Renderer::display_image(const std::string& filename) {
    stbi_write_png(filename.c_str(), W, H, 3, image.data(), 0);
}


void Renderer::render(const std::string& filename) {
    const int total_pixels = W * H;
    static std::atomic<int> pixel_count{0};
    pixel_count = 0; 

    auto start_time = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Ray ray = camera.get_ray(i, j);
            Intersection intersection = scene.get_closest_hit(ray);
            if (intersection.flag) {
                Vector vec_albedo_average;
                for (int k = 0; k < nr_rays; k++) {
                    Vector vec_albedo = scene.get_intensity(ray, intensity_depth, false);
                    vec_albedo_average = vec_albedo_average + vec_albedo;
                }
                vec_albedo_average = (1.0 / (static_cast<double>(nr_rays))) * vec_albedo_average;
                set_image_pixel_color(i, j, vec_albedo_average);
            }

            int local_count = ++pixel_count;
            if (local_count % 5000 == 0) {
                double fraction = double(local_count) / double(total_pixels);

                const int barWidth = 50;
                int fill = static_cast<int>(fraction * barWidth);

                #pragma omp critical
                {
                    std::cout << "\r[";
                    for (int b = 0; b < barWidth; ++b) {
                        if (b < fill)
                            std::cout << "=";  
                        else
                            std::cout << " ";  
                    }
                    std::cout << "] " 
                              << int(fraction * 100.0) << "%  " 
                              << std::flush;  
                }
            }
        }
    }

    display_image(filename);

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    std::cout << std::endl << "Render loop took " << elapsed.count() << " seconds.\n";
}
