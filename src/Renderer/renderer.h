#ifndef RENDERER_H
#define RENDERER_H

#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <atomic>

#include "vector.h"
#include "ray.h"
#include "sphere.h"
#include "scene.h"
#include "camera.h"

#define INTENSITY_DEPTH 8

class Renderer {
    public:
        explicit Renderer(
            Scene aux_scene = Scene(), 
            Camera aux_camera = Camera(),
            int aux_nr_rays = 1
        );

        void set_image_pixel_color(const int& i, const int& j, Vector& vec_albedo);
        Vector get_image_pixel_color(const int& i, const int& j);

        void display_image();

        void render();

    private:
        Scene scene;
        Camera camera;
        int W;
        int H;
        int nr_rays;

        std::vector<unsigned char> image; // The image that will be rendered
    };

#endif // RENDERER_H
