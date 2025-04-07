#ifndef RENDERER_H
#define RENDERER_H

#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <algorithm>
#include <cmath>

#include "vector.h"
#include "ray.h"
#include "sphere.h"
#include "scene.h"
#include "camera.h"


class Renderer {
    public:
        explicit Renderer(
            Scene aux_scene = Scene(), 
            Camera aux_camera = Camera()
        );

        void set_image_pixel_color(int i, int j, Vector vec_albedo);
        Vector get_image_pixel_color(int i, int j);

        void display_image();

        void render_default_image();
        void render_ray_scene_intersection();
        void render_shadows();
        void render_reflections_shadows();
        void render_refractions_reflections_shadows();

    private:
        Scene scene;
        Camera camera;
        int W;
        int H;

        std::vector<unsigned char> image; // The image that will be rendered
    };

#endif // RENDERER_H
