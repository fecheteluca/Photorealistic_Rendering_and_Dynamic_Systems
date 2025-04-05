#include "scene.h"

// Scene class definitions
Scene::Scene(std::vector<Sphere> aux_l_sph, Vector aux_vec_light_source, double aux_light_intensity) { 
    l_sph = aux_l_sph;
    vec_light_source = aux_vec_light_source;
    light_intensity = aux_light_intensity;
}

void Scene::add_object(Sphere sph_extra) {
    l_sph.push_back(sph_extra);
}

Intersection Scene::get_closest_hit(Ray ray) {
    double closest_hit_distance = INT_MAX * 1.0;
    Intersection closest_hit_intersection = Intersection();
    
    for (auto sph = l_sph.begin(); sph != l_sph.end(); sph++) {
        Intersection intersection = (*sph).intersected_by(ray);
        if (intersection.flag && intersection.distance < closest_hit_distance) {
            closest_hit_distance = intersection.distance;
            closest_hit_intersection = intersection;
        }
    }

    return closest_hit_intersection;
}

Vector Scene::get_shadow_intensity(Intersection intersection) {
    double distance = (vec_light_source - intersection.vec_point).norm();
    Vector omega    = vec_light_source - intersection.vec_point;
    omega.normalize();

    double first_factor = light_intensity / (4 * M_PI * pow(distance, 2));
    Vector second_factor = (1 / M_PI) * intersection.vec_albedo;

    double third_factor = 0;
    Vector vec_offseted_point = intersection.vec_point + 0.001 * intersection.vec_normal;
    Ray ray = Ray(vec_offseted_point, omega);   
    Intersection src_pnt_intersection = get_closest_hit(ray);
    if (distance <= src_pnt_intersection.distance) {
        third_factor = 1;
    }

    double fourth_factor = dot(intersection.vec_normal, omega);
    return first_factor * second_factor * third_factor * fourth_factor;
}


// The definition of the standard scene that we are going to use throughout the entire project.
Scene get_standard_scene() {
    std::vector<Sphere> l_sph_sc_standard(7);

    // Center Sphere Configuration
    Vector sph_center_vec_center = Vector(0, 0, 0);
    double sph_center_radius     = 10;
    Vector sph_center_vec_albedo = Vector(0.8, 0.8, 0.8);
    Sphere sph_center = Sphere(sph_center_vec_center, sph_center_radius, sph_center_vec_albedo);
    l_sph_sc_standard.push_back(sph_center);

    // Left Wall Sphere Configuration
    Vector sph_leftwall_vec_center = Vector(1000, 0, 0);
    double sph_leftwall_radius     = 940;
    Vector sph_leftwall_vec_albedo = Vector(0.9, 0.2, 0.9);
    Sphere sph_leftwall = Sphere(sph_leftwall_vec_center, sph_leftwall_radius, sph_leftwall_vec_albedo);
    l_sph_sc_standard.push_back(sph_leftwall);

    // Right Wall Sphere Configuration
    Vector sph_rightwall_vec_center = Vector(-1000, 0, 0);
    double sph_rightwall_radius     = 940;
    Vector sph_rightwall_vec_albedo = Vector(0.6, 0.5, 0.1);
    Sphere sph_rightwall = Sphere(sph_rightwall_vec_center, sph_rightwall_radius, sph_rightwall_vec_albedo);
    l_sph_sc_standard.push_back(sph_rightwall);

    // Up Wall (Ceiling) Sphere Configuration
    Vector sph_upwall_vec_center = Vector(0, 1000, 0);
    double sph_upwall_radius     = 940;
    Vector sph_upwall_vec_albedo = Vector(0.2, 0.5, 0.9);
    Sphere sph_upwall = Sphere(sph_upwall_vec_center, sph_upwall_radius, sph_upwall_vec_albedo);
    l_sph_sc_standard.push_back(sph_upwall);

    // Down Wall (Floor) Sphere Configuration
    Vector sph_downwall_vec_center = Vector(0, -1000, 0);
    double sph_downwall_radius     = 990;
    Vector sph_downwall_vec_albedo = Vector(0.3, 0.4, 0.7);
    Sphere sph_downwall = Sphere(sph_downwall_vec_center, sph_downwall_radius, sph_downwall_vec_albedo);
    l_sph_sc_standard.push_back(sph_downwall);

    // Back Wall Sphere Configuration
    Vector sph_backwall_vec_center = Vector(0, 0, 1000);
    double sph_backwall_radius     = 940;
    Vector sph_backwall_vec_albedo = Vector(0.9, 0.4, 0.3);
    Sphere sph_backwall = Sphere(sph_backwall_vec_center, sph_backwall_radius, sph_backwall_vec_albedo);
    l_sph_sc_standard.push_back(sph_backwall);

    // Front Wall Sphere Configuration
    Vector sph_frontwall_vec_center = Vector(0, 0, -1000);
    double sph_frontwall_radius     = 940;
    Vector sph_frontwall_vec_albedo = Vector(0.4, 0.8, 0.7);
    Sphere sph_frontwall = Sphere(sph_frontwall_vec_center, sph_frontwall_radius, sph_frontwall_vec_albedo);
    l_sph_sc_standard.push_back(sph_frontwall);

    Vector vec_light_source = Vector(-10, 20, 40);
    double light_intensity = 40000;

    return Scene(l_sph_sc_standard, vec_light_source, light_intensity);
}