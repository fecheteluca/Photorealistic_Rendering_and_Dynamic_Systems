#include "scene.h"

// Scene class definitions
Scene::Scene(
    std::vector<Sphere> aux_l_sph, 
    Vector aux_vec_light_source, 
    double aux_light_intensity,
    double aux_refraction_index
) { 
    l_sph = aux_l_sph;
    vec_light_source = aux_vec_light_source;
    light_intensity = aux_light_intensity;
    refraction_index = aux_refraction_index;
}

void Scene::add_object(Sphere sph_extra) {
    l_sph.push_back(sph_extra);
}

Intersection Scene::get_closest_hit(Ray ray) {
    double closest_hit_distance = INT_MAX * 1.0;
    Intersection closest_hit_intersection = Intersection();
    
    for (int idx = 0; idx <= l_sph.size(); idx++) {
        Intersection intersection = l_sph[idx].intersected_by(ray);
        if (intersection.flag && intersection.distance < closest_hit_distance) {
            closest_hit_distance = intersection.distance;
            closest_hit_intersection = intersection;
            closest_hit_intersection.idx_sph = idx;
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

Vector Scene::get_reflected_intensity(Ray ray , int ray_depth) {
    if (ray_depth < 0) {
        return Vector(0.0, 0.0, 0.0);
    }
    else {
        Intersection first_hit_intersection = get_closest_hit(ray);
        if (first_hit_intersection.flag) {
            if (l_sph[first_hit_intersection.idx_sph].has_mirror_surface()) {
                Vector vec_unit_direction = ray.get_unit_direction();
                Vector vec_normal_surface = first_hit_intersection.vec_normal;

                Vector vec_reflected_origin = first_hit_intersection.vec_point + 0.001 * vec_normal_surface;

                Vector vec_reflected_unit_direction =  vec_unit_direction - 2 * dot(vec_unit_direction, vec_normal_surface) * vec_normal_surface;
                vec_reflected_unit_direction.normalize();

                Ray ray_reflected = Ray(vec_reflected_origin, vec_reflected_unit_direction);

                return get_reflected_intensity(ray_reflected, ray_depth - 1);
            }
            else {
                return get_shadow_intensity(first_hit_intersection);
            }
        }
    }
}

Vector Scene::get_refracted_intensity(Ray ray , int ray_depth) {
    if (ray_depth < 0) {
        return Vector(0.0, 0.0, 0.0);
    }
    else {
        Intersection first_hit_intersection = get_closest_hit(ray);
        if (first_hit_intersection.flag) {
            if (l_sph[first_hit_intersection.idx_sph].has_mirror_surface()) {
                Vector vec_unit_direction = ray.get_unit_direction();
                Vector vec_normal_surface = first_hit_intersection.vec_normal;

                Vector vec_reflected_origin = first_hit_intersection.vec_point + 0.001 * vec_normal_surface;

                Vector vec_reflected_unit_direction =  vec_unit_direction - 2 * dot(vec_unit_direction, vec_normal_surface) * vec_normal_surface;
                vec_reflected_unit_direction.normalize();

                Ray ray_reflected = Ray(vec_reflected_origin, vec_reflected_unit_direction);

                return get_refracted_intensity(ray_reflected, ray_depth - 1);
            }
            else if (l_sph[first_hit_intersection.idx_sph].has_transparent_surface()){
                Vector vec_unit_direction = ray.get_unit_direction();
                Vector vec_normal_surface = first_hit_intersection.vec_normal;
                double cosine_unit_normal = dot(vec_unit_direction, vec_normal_surface);

                double sphere_refraction_index = l_sph[first_hit_intersection.idx_sph].get_refraction_index();
                double scene_refraction_index = refraction_index;

                double eta_1 = cosine_unit_normal <= 0 ? scene_refraction_index : sphere_refraction_index;
                double eta_2 = cosine_unit_normal <= 0 ? sphere_refraction_index : scene_refraction_index;
                double refraction_ratio = eta_1 / eta_2;

                vec_normal_surface = cosine_unit_normal <= 0 ? vec_normal_surface : (-1) * vec_normal_surface;
                cosine_unit_normal = cosine_unit_normal <= 0 ? cosine_unit_normal : (-1) * cosine_unit_normal;

                double k = 1 - pow(refraction_ratio, 2) * (1 - pow(cosine_unit_normal, 2));
            
                if (k < 0) {
                    Vector vec_reflected_origin = first_hit_intersection.vec_point + 0.001 * vec_normal_surface;

                    Vector vec_reflected_unit_direction = vec_unit_direction - 2 * dot(vec_unit_direction, vec_normal_surface) * vec_normal_surface;
                    vec_reflected_unit_direction.normalize();

                    Ray ray_reflected = Ray(vec_reflected_origin, vec_reflected_unit_direction);

                    return get_refracted_intensity(ray_reflected, ray_depth - 1);
                }
                
                Vector vec_refracted_origin = first_hit_intersection.vec_point - 0.001 * vec_normal_surface;

                Vector vec_refracted_unit_direction_tangential = refraction_ratio * (vec_unit_direction - cosine_unit_normal * vec_normal_surface);
                Vector vec_refracted_unit_direction_normal = (-1) * vec_normal_surface * sqrt(k);
                Vector vec_refracted_unit_direction = vec_refracted_unit_direction_tangential + vec_refracted_unit_direction_normal;
                vec_refracted_unit_direction.normalize();
                    
                Ray ray_refracted(vec_refracted_origin, vec_refracted_unit_direction);

                return get_refracted_intensity(ray_refracted, ray_depth - 1);
            }
            else {
                return get_shadow_intensity(first_hit_intersection);
            }
        }
    }
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

Scene get_intermediate_scene_reflection() {
    std::vector<Sphere> l_sph_sc_intermediate(9);

    // Center Sphere Configuration
    Vector sph_center_vec_center = Vector(0, 0, 0);
    double sph_center_radius     = 10;
    Vector sph_center_vec_albedo = Vector(0.8, 0.8, 0.8);
    bool sph_center_mirror = true;
    Sphere sph_center = Sphere(sph_center_vec_center, sph_center_radius, sph_center_vec_albedo, sph_center_mirror);
    l_sph_sc_intermediate.push_back(sph_center);

    // Center Left Sphere Configuration
    Vector sph_center_left_vec_center = Vector(20, 0, 0);
    double sph_center_left_radius     = 10;
    Vector sph_center_left_vec_albedo = Vector(0.8, 0.8, 0.8);
    bool sph_center_left_mirror = true;
    Sphere sph_center_left = Sphere(sph_center_left_vec_center, sph_center_left_radius, sph_center_left_vec_albedo, sph_center_left_mirror);
    l_sph_sc_intermediate.push_back(sph_center_left);

    // Center Right Sphere Configuration
    Vector sph_center_right_vec_center = Vector(-20, 0, 0);
    double sph_center_right_radius     = 10;
    Vector sph_center_right_vec_albedo = Vector(0.8, 0.8, 0.8);
    bool sph_center_right_mirror = true;
    Sphere sph_center_right = Sphere(sph_center_right_vec_center, sph_center_right_radius, sph_center_right_vec_albedo, sph_center_right_mirror);
    l_sph_sc_intermediate.push_back(sph_center_right);

    // Left Wall Sphere Configuration
    Vector sph_leftwall_vec_center = Vector(1000, 0, 0);
    double sph_leftwall_radius     = 940;
    Vector sph_leftwall_vec_albedo = Vector(0.9, 0.2, 0.9);
    Sphere sph_leftwall = Sphere(sph_leftwall_vec_center, sph_leftwall_radius, sph_leftwall_vec_albedo);
    l_sph_sc_intermediate.push_back(sph_leftwall);

    // Right Wall Sphere Configuration
    Vector sph_rightwall_vec_center = Vector(-1000, 0, 0);
    double sph_rightwall_radius     = 940;
    Vector sph_rightwall_vec_albedo = Vector(0.6, 0.5, 0.1);
    Sphere sph_rightwall = Sphere(sph_rightwall_vec_center, sph_rightwall_radius, sph_rightwall_vec_albedo);
    l_sph_sc_intermediate.push_back(sph_rightwall);

    // Up Wall (Ceiling) Sphere Configuration
    Vector sph_upwall_vec_center = Vector(0, 1000, 0);
    double sph_upwall_radius     = 940;
    Vector sph_upwall_vec_albedo = Vector(0.2, 0.5, 0.9);
    Sphere sph_upwall = Sphere(sph_upwall_vec_center, sph_upwall_radius, sph_upwall_vec_albedo);
    l_sph_sc_intermediate.push_back(sph_upwall);

    // Down Wall (Floor) Sphere Configuration
    Vector sph_downwall_vec_center = Vector(0, -1000, 0);
    double sph_downwall_radius     = 990;
    Vector sph_downwall_vec_albedo = Vector(0.3, 0.4, 0.7);
    Sphere sph_downwall = Sphere(sph_downwall_vec_center, sph_downwall_radius, sph_downwall_vec_albedo);
    l_sph_sc_intermediate.push_back(sph_downwall);

    // Back Wall Sphere Configuration
    Vector sph_backwall_vec_center = Vector(0, 0, 1000);
    double sph_backwall_radius     = 940;
    Vector sph_backwall_vec_albedo = Vector(0.9, 0.4, 0.3);
    Sphere sph_backwall = Sphere(sph_backwall_vec_center, sph_backwall_radius, sph_backwall_vec_albedo);
    l_sph_sc_intermediate.push_back(sph_backwall);

    // Front Wall Sphere Configuration
    Vector sph_frontwall_vec_center = Vector(0, 0, -1000);
    double sph_frontwall_radius     = 940;
    Vector sph_frontwall_vec_albedo = Vector(0.4, 0.8, 0.7);
    Sphere sph_frontwall = Sphere(sph_frontwall_vec_center, sph_frontwall_radius, sph_frontwall_vec_albedo);
    l_sph_sc_intermediate.push_back(sph_frontwall);

    Vector vec_light_source = Vector(-10, 20, 40);
    double light_intensity = 40000;

    return Scene(l_sph_sc_intermediate, vec_light_source, light_intensity);
}

Scene get_intermediate_scene_refraction_reflection() {
    std::vector<Sphere> l_sph_sc_intermediate(9);

    // Center Sphere Configuration
    Vector sph_center_vec_center = Vector(0, 0, 0);
    double sph_center_radius     = 10;
    Vector sph_center_vec_albedo = Vector(0.8, 0.8, 0.8);
    bool sph_center_mirror = false;
    bool sph_center_transparent = true;
    double sph_center_refraction_index = 1.5;
    Sphere sph_center = Sphere(sph_center_vec_center, sph_center_radius, sph_center_vec_albedo, sph_center_mirror, sph_center_transparent, sph_center_refraction_index);
    l_sph_sc_intermediate.push_back(sph_center);

    // Center Left Sphere Configuration
    Vector sph_center_left_vec_center = Vector(-20, 0, 0);
    double sph_center_left_radius     = 10;
    Vector sph_center_left_vec_albedo = Vector(0.8, 0.8, 0.8);
    bool sph_center_left_mirror = true;
    bool sph_center_left_transparent = false;
    Sphere sph_center_left = Sphere(sph_center_left_vec_center, sph_center_left_radius, sph_center_left_vec_albedo, sph_center_left_mirror, sph_center_left_transparent);
    l_sph_sc_intermediate.push_back(sph_center_left);

    // Center Right Sphere Configuration
    Vector sph_center_right_vec_center = Vector(20, 0, 0);
    double sph_center_right_radius     = 10;
    Vector sph_center_right_vec_albedo = Vector(0.8, 0.8, 0.8);
    bool sph_center_right_mirror = true;
    bool sph_center_right_transparent = false;
    Sphere sph_center_right = Sphere(sph_center_right_vec_center, sph_center_right_radius, sph_center_right_vec_albedo, sph_center_right_mirror, sph_center_right_transparent);
    l_sph_sc_intermediate.push_back(sph_center_right);

    // Left Wall Sphere Configuration
    Vector sph_leftwall_vec_center = Vector(1000, 0, 0);
    double sph_leftwall_radius     = 940;
    Vector sph_leftwall_vec_albedo = Vector(0.9, 0.2, 0.9);
    Sphere sph_leftwall = Sphere(sph_leftwall_vec_center, sph_leftwall_radius, sph_leftwall_vec_albedo);
    l_sph_sc_intermediate.push_back(sph_leftwall);

    // Right Wall Sphere Configuration
    Vector sph_rightwall_vec_center = Vector(-1000, 0, 0);
    double sph_rightwall_radius     = 940;
    Vector sph_rightwall_vec_albedo = Vector(0.6, 0.5, 0.1);
    Sphere sph_rightwall = Sphere(sph_rightwall_vec_center, sph_rightwall_radius, sph_rightwall_vec_albedo);
    l_sph_sc_intermediate.push_back(sph_rightwall);

    // Up Wall (Ceiling) Sphere Configuration
    Vector sph_upwall_vec_center = Vector(0, 1000, 0);
    double sph_upwall_radius     = 940;
    Vector sph_upwall_vec_albedo = Vector(0.2, 0.5, 0.9);
    Sphere sph_upwall = Sphere(sph_upwall_vec_center, sph_upwall_radius, sph_upwall_vec_albedo);
    l_sph_sc_intermediate.push_back(sph_upwall);

    // Down Wall (Floor) Sphere Configuration
    Vector sph_downwall_vec_center = Vector(0, -1000, 0);
    double sph_downwall_radius     = 990;
    Vector sph_downwall_vec_albedo = Vector(0.3, 0.4, 0.7);
    Sphere sph_downwall = Sphere(sph_downwall_vec_center, sph_downwall_radius, sph_downwall_vec_albedo);
    l_sph_sc_intermediate.push_back(sph_downwall);

    // Back Wall Sphere Configuration
    Vector sph_backwall_vec_center = Vector(0, 0, 1000);
    double sph_backwall_radius     = 940;
    Vector sph_backwall_vec_albedo = Vector(0.9, 0.4, 0.3);
    Sphere sph_backwall = Sphere(sph_backwall_vec_center, sph_backwall_radius, sph_backwall_vec_albedo);
    l_sph_sc_intermediate.push_back(sph_backwall);

    // Front Wall Sphere Configuration
    Vector sph_frontwall_vec_center = Vector(0, 0, -1000);
    double sph_frontwall_radius     = 940;
    Vector sph_frontwall_vec_albedo = Vector(0.4, 0.8, 0.7);
    Sphere sph_frontwall = Sphere(sph_frontwall_vec_center, sph_frontwall_radius, sph_frontwall_vec_albedo);
    l_sph_sc_intermediate.push_back(sph_frontwall);

    Vector vec_light_source = Vector(-10, 20, 40);
    double light_intensity = 40000;
    double scene_refraction_index = 1.0;

    return Scene(l_sph_sc_intermediate, vec_light_source, light_intensity, scene_refraction_index);
}