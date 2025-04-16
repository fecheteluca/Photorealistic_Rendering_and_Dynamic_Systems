#include "renderer.h"

int main() {
    std::vector<Sphere> l_sph(11);

    // Center Sphere Configuration
    Vector sph_center_vec_center = Vector(0, 0, 0);
    double sph_center_radius     = 10;
    Vector sph_center_vec_albedo = Vector(0.8, 0.8, 0.8);
    bool sph_center_mirror = false;
    bool sph_center_transparent = true;
    bool sph_center_light_source = false;
    double sph_center_refraction_index = 1.5;
    bool sph_center_invert_normals = false;
    Sphere sph_center = Sphere(
        sph_center_vec_center, 
        sph_center_radius, 
        sph_center_vec_albedo, 
        sph_center_mirror, 
        sph_center_transparent, 
        sph_center_light_source,
        sph_center_refraction_index,
        sph_center_invert_normals
    );
    l_sph.push_back(sph_center);

    // Center Left Sphere Configuration
    Vector sph_center_left_vec_center = Vector(-20, 0, 0);
    double sph_center_left_radius     = 10;
    Vector sph_center_left_vec_albedo = Vector(0.8, 0.8, 0.8);
    bool sph_center_left_mirror = true;
    bool sph_center_left_transparent = false;
    bool sph_center_left_light_source = false;
    double sph_center_left_refraction_index = 1.0;
    bool sph_center_left_invert_normals = false;
    Sphere sph_center_left = Sphere(
        sph_center_left_vec_center, 
        sph_center_left_radius, 
        sph_center_left_vec_albedo, 
        sph_center_left_mirror, 
        sph_center_left_transparent,
        sph_center_left_light_source,
        sph_center_left_refraction_index,
        sph_center_left_invert_normals
    );
    l_sph.push_back(sph_center_left);

    // Center Right Outside Sphere Configuration
    Vector sph_center_right_outside_vec_center = Vector(20, 0, 0);
    double sph_center_right_outside_radius     = 10;
    Vector sph_center_right_outside_vec_albedo = Vector(0.8, 0.8, 0.8);
    bool sph_center_right_outside_mirror = false;
    bool sph_center_right_outside_transparent = true;
    bool sph_center_right_outside_light_source = false;
    double sph_center_right_outside_refraction_index = 1.5;
    bool sph_center_right_outside_invert_normals = false;
    Sphere sph_center_right_outside = Sphere(
        sph_center_right_outside_vec_center, 
        sph_center_right_outside_radius, 
        sph_center_right_outside_vec_albedo, 
        sph_center_right_outside_mirror, 
        sph_center_right_outside_transparent,
        sph_center_right_outside_light_source,
        sph_center_right_outside_refraction_index,
        sph_center_right_outside_invert_normals
    );
    l_sph.push_back(sph_center_right_outside);

    // Center Right Inside Sphere Configuration
    Vector sph_center_right_inside_vec_center = Vector(20, 0, 0);
    double sph_center_right_inside_radius     = 9.6;
    Vector sph_center_right_inside_vec_albedo = Vector(0.8, 0.8, 0.8);
    bool sph_center_right_inside_mirror = false;
    bool sph_center_right_inside_transparent = true;
    bool sph_center_right_inside_light_source = false;
    double sph_center_right_inside_refraction_index = 1.5;
    bool sph_center_right_inside_invert_normals = true;
    Sphere sph_center_right_inside = Sphere(
        sph_center_right_inside_vec_center, 
        sph_center_right_inside_radius, 
        sph_center_right_inside_vec_albedo, 
        sph_center_right_inside_mirror, 
        sph_center_right_inside_transparent,
        sph_center_right_inside_light_source,
        sph_center_right_inside_refraction_index,
        sph_center_right_inside_invert_normals
    );
    l_sph.push_back(sph_center_right_inside);

    // Left Wall Sphere Configuration
    Vector sph_leftwall_vec_center = Vector(1000, 0, 0);
    double sph_leftwall_radius     = 940;
    Vector sph_leftwall_vec_albedo = Vector(0.9, 0.2, 0.9);
    Sphere sph_leftwall = Sphere(
        sph_leftwall_vec_center, 
        sph_leftwall_radius, 
        sph_leftwall_vec_albedo
    );
    l_sph.push_back(sph_leftwall);

    // Right Wall Sphere Configuration
    Vector sph_rightwall_vec_center = Vector(-1000, 0, 0);
    double sph_rightwall_radius     = 940;
    Vector sph_rightwall_vec_albedo = Vector(0.6, 0.5, 0.1);
    Sphere sph_rightwall = Sphere(
        sph_rightwall_vec_center, 
        sph_rightwall_radius, 
        sph_rightwall_vec_albedo
    );
    l_sph.push_back(sph_rightwall);

    // Up Wall (Ceiling) Sphere Configuration
    Vector sph_upwall_vec_center = Vector(0, 1000, 0);
    double sph_upwall_radius     = 940;
    Vector sph_upwall_vec_albedo = Vector(0.2, 0.5, 0.9);
    Sphere sph_upwall = Sphere(
        sph_upwall_vec_center, 
        sph_upwall_radius, 
        sph_upwall_vec_albedo
    );
    l_sph.push_back(sph_upwall);

    // Down Wall (Floor) Sphere Configuration
    Vector sph_downwall_vec_center = Vector(0, -1000, 0);
    double sph_downwall_radius     = 990;
    Vector sph_downwall_vec_albedo = Vector(0.3, 0.4, 0.7);
    Sphere sph_downwall = Sphere(
        sph_downwall_vec_center, 
        sph_downwall_radius, 
        sph_downwall_vec_albedo
    );
    l_sph.push_back(sph_downwall);

    // Back Wall Sphere Configuration
    Vector sph_backwall_vec_center = Vector(0, 0, 1000);
    double sph_backwall_radius     = 940;
    Vector sph_backwall_vec_albedo = Vector(0.9, 0.4, 0.3);
    Sphere sph_backwall = Sphere(
        sph_backwall_vec_center, 
        sph_backwall_radius, 
        sph_backwall_vec_albedo
    );
    l_sph.push_back(sph_backwall);

    // Front Wall Sphere Configuration
    Vector sph_frontwall_vec_center = Vector(0, 0, -1000);
    double sph_frontwall_radius     = 940;
    Vector sph_frontwall_vec_albedo = Vector(0.4, 0.8, 0.7);
    Sphere sph_frontwall = Sphere(
        sph_frontwall_vec_center, 
        sph_frontwall_radius, 
        sph_frontwall_vec_albedo
    );
    l_sph.push_back(sph_frontwall);

    // Light Source Sphere Configuration
    Vector sph_light_source_vec_center = Vector(-10, 25, -10);
    double sph_light_source_radius     = 5;
    Vector sph_light_source_vec_albedo = Vector(0.0, 0.0, 0.0);
    bool sph_light_source_mirror = false;
    bool sph_light_source_transparent = false;
    bool sph_light_source_light_source = true;
    double sph_light_source_refraction_index = 1.5;
    bool sph_light_source_invert_normals = false;
    Sphere sph_light_source = Sphere(
        sph_light_source_vec_center, 
        sph_light_source_radius, 
        sph_light_source_vec_albedo, 
        sph_light_source_mirror, 
        sph_light_source_transparent, 
        sph_light_source_light_source,
        sph_light_source_refraction_index,
        sph_light_source_invert_normals
    );

    bool light_source_sphere = false;

    if (light_source_sphere) {
        l_sph.push_back(sph_light_source);
    }
    
    Vector vec_light_source = Vector(-10, 20, 40);
    double light_intensity = 1e05;
    double scene_refraction_index = 1.0;

    Scene scene = Scene(
        l_sph, 
        vec_light_source, 
        sph_light_source,
        light_intensity, 
        scene_refraction_index,
        light_source_sphere
    );

    Vector vec_camera_center = Vector(0, 0, 55);
    double camera_angle = 60;
    int camera_width = 512;
    int camera_height = 512;
    double camera_stdev = 1.0;
    double camera_spread = 0.5;
    double camera_focal_dist = 55;
    double camera_aperture = 1.0;
    Camera camera = Camera(
        vec_camera_center, 
        camera_angle, 
        camera_width, 
        camera_height, 
        camera_stdev, 
        camera_spread, 
        camera_focal_dist, 
        camera_aperture
    );

    int rays_per_pixel = 500;
    int intensity_depth = 8;
    Renderer renderer = Renderer(
        scene, 
        camera, 
        rays_per_pixel, 
        intensity_depth
    );
    renderer.render("../../images/image.png");

    return 0;
}
