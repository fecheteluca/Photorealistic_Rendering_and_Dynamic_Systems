#ifndef RENDERER_H
#define RENDERER_H

#include <vector>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <atomic>
#include <sstream>

#include "vector.h"
#include "ray.h"
#include "sphere.h"
#include "scene.h"
#include "camera.h"
#include "triangle_mesh.h"
#include "stb_image_write.h"
#include "polygon.h"
#include "optimal_transport.h"
#include "particle_simulation.h"

class Renderer {
    public:
        explicit Renderer(
            Scene aux_scene = Scene(), 
            Camera aux_camera = Camera(),
            int aux_nr_rays = 1,
            int aux_intensity_depth = 5
        );

        void set_image_pixel_color(const int& i, const int& j, Vector& vec_albedo);
        Vector get_image_pixel_color(const int& i, const int& j);

        void display_image(const std::string& filename);

        void render(const std::string& filename);

    private:
        Scene scene;
        Camera camera;
        int W;
        int H;
        int nr_rays;
        int intensity_depth;

        std::vector<unsigned char> image; // The image that will be rendered
    };

    // saves a static svg file. The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
    void save_svg(const std::vector<Polygon> &polygons, std::string filename, std::string fillcol = "none");
 
 
    // Adds one frame of an animated svg file. frameid is the frame number (between 0 and nbframes-1).
    // polygons is a list of polygons, describing the current frame.
    // The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
    void save_svg_animated(const std::vector<Polygon> &polygons,
                       const std::vector<std::string>& colors,
                       std::string filename, int frameid, int nbframes);

    void test_polygon_clipping();
    void test_and_visualize_voronoi();
    void run_and_visualize_power_diagram();
    void run_and_optimize_power_diagram();
    void run_and_visualize_particle_simulation();

#endif // RENDERER_H
