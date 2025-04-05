#ifndef RENDERER_H
#define RENDERER_H

#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <cmath>

#include "vector.h"
#include "ray.h"
#include "sphere.h"
#include "scene.h"
#include "camera.h"


class Renderer {
    public:
        explicit Renderer(Scene aux_scene = Scene(), Camera aux_camera = Camera(), Vector aux_vec_light_source = Vector(), double light_intensity = 0);

        void set_image_pixel_color(int i, int j, Vector vec_albedo);
        Vector get_image_pixel_color(int i, int j);
        Vector get_reflected_intensity(Intersection intersection);

        void display_image();

        void render_default_image();
        void render_ray_scene_intersection();
        void render_shading_and_shadows();

    private:
        Scene scene;
        Camera camera;
        int W;
        int H;
        Vector vec_light_source;
        double light_intensity;

        std::vector<unsigned char> image; // The image that will be rendered
    };

#endif // RENDERER_H
