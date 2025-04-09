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

void Scene::add_object(const Sphere& sph_extra) {
    l_sph.push_back(sph_extra);
}

Intersection Scene::get_closest_hit(Ray& ray) {
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

Vector Scene::get_shadow_intensity(const Intersection& intersection) {
    double distance = (vec_light_source - intersection.vec_point).norm();
    Vector omega    = vec_light_source - intersection.vec_point;
    omega.normalize();

    double first_factor = light_intensity / (4 * M_PI * pow(distance, 2));
    Vector second_factor = (1 / M_PI) * intersection.vec_albedo;

    double third_factor = 0;
    Vector vec_offseted_point = intersection.vec_point + EPS * intersection.vec_normal;
    Ray ray = Ray(vec_offseted_point, omega);   
    Intersection src_pnt_intersection = get_closest_hit(ray);
    if (distance <= src_pnt_intersection.distance) {
        third_factor = 1;
    }

    double fourth_factor = dot(intersection.vec_normal, omega);
    return first_factor * second_factor * third_factor * fourth_factor;
}  

Ray Scene::random_cos(const Intersection& intersection) {
    double r_1 = uniform(engine);
    double r_2 = uniform(engine);

    double x = cos(2 * M_PI * r_1) * sqrt(1 - r_2);
    double y = sin(2 * M_PI * r_1) * sqrt(1 - r_2);;
    double z = sqrt(r_2);

    Vector vec_normal = intersection.vec_normal;
    double x_vec_normal = vec_normal.get_x();
    double y_vec_normal = vec_normal.get_y();
    double z_vec_normal = vec_normal.get_z();

    double smallest_comp = std::min(
        fabs(x_vec_normal),
        std::min(fabs(y_vec_normal), fabs(z_vec_normal))
    );    

    Vector T_1 = Vector();
    if (smallest_comp == fabs(x_vec_normal)) {
        T_1 = Vector(0.0, z_vec_normal, -y_vec_normal);
    }
    else if (smallest_comp == fabs(y_vec_normal)) {
        T_1 = Vector(-z_vec_normal, 0.0, x_vec_normal);
    }
    else if (smallest_comp == fabs(z_vec_normal)) {
        T_1 = Vector(y_vec_normal, -x_vec_normal, 0.0);
    }
    T_1.normalize();

    Vector T_2 = cross(vec_normal, T_1);

    Vector vec_random_origin = intersection.vec_point + EPS * vec_normal;
    Vector vec_random_unit_direction = x * T_1 + y * T_2 + z * vec_normal;
    return Ray(vec_random_origin, vec_random_unit_direction);
}

Vector Scene::get_intensity(Ray& ray , const int& ray_depth) {
    if (ray_depth < 0) {
        return Vector(0.0, 0.0, 0.0);
    }
    else {
        Intersection first_hit_intersection = get_closest_hit(ray);
        if (first_hit_intersection.flag) {
            if (l_sph[first_hit_intersection.idx_sph].has_mirror_surface()) {
                Vector vec_unit_direction = ray.get_unit_direction();
                Vector vec_normal_surface = first_hit_intersection.vec_normal;

                Vector vec_reflected_origin = first_hit_intersection.vec_point + EPS * vec_normal_surface;

                Vector vec_reflected_unit_direction =  vec_unit_direction - 2 * dot(vec_unit_direction, vec_normal_surface) * vec_normal_surface;
                vec_reflected_unit_direction.normalize();

                Ray ray_reflected = Ray(vec_reflected_origin, vec_reflected_unit_direction);

                return get_intensity(ray_reflected, ray_depth - 1);
            }
            else if (l_sph[first_hit_intersection.idx_sph].has_transparent_surface()){
                Vector vec_unit_direction = ray.get_unit_direction();
                Vector vec_normal_surface = first_hit_intersection.vec_normal;
                
                // Reflected ray computation to be used for the Fresnel Law
                Vector vec_reflected_origin = first_hit_intersection.vec_point + EPS * vec_normal_surface;

                Vector vec_reflected_unit_direction = vec_unit_direction - 2 * dot(vec_unit_direction, vec_normal_surface) * vec_normal_surface;
                vec_reflected_unit_direction.normalize();

                Ray ray_reflected = Ray(vec_reflected_origin, vec_reflected_unit_direction);

                // Refracted ray computation to be used for the Fresnel Law
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
                    Vector vec_reflected_origin = first_hit_intersection.vec_point + EPS * vec_normal_surface;

                    Vector vec_reflected_unit_direction = vec_unit_direction - 2 * dot(vec_unit_direction, vec_normal_surface) * vec_normal_surface;
                    vec_reflected_unit_direction.normalize();

                    Ray ray_reflected = Ray(vec_reflected_origin, vec_reflected_unit_direction);

                    return get_intensity(ray_reflected, ray_depth - 1);
                }
                
                Vector vec_refracted_origin = first_hit_intersection.vec_point - EPS * vec_normal_surface;

                Vector vec_refracted_unit_direction_tangential = refraction_ratio * (vec_unit_direction - cosine_unit_normal * vec_normal_surface);
                Vector vec_refracted_unit_direction_normal = (-1) * vec_normal_surface * sqrt(k);
                Vector vec_refracted_unit_direction = vec_refracted_unit_direction_tangential + vec_refracted_unit_direction_normal;
                vec_refracted_unit_direction.normalize();
                    
                Ray ray_refracted = Ray(vec_refracted_origin, vec_refracted_unit_direction);

                // Fresnel Law Equations
                double k_0 = pow((eta_1 - eta_2), 2) / pow((eta_1 + eta_2), 2);
                double R   = k_0 + (1 - k_0) * pow((1 - fabs(cosine_unit_normal)), 5);

                double u = uniform(engine);
                Ray ray_returned = u < R ? ray_reflected : ray_refracted;

                return get_intensity(ray_returned, ray_depth - 1);
            }
            else {
                Vector Lo = get_shadow_intensity(first_hit_intersection);
                Ray random_ray = random_cos(first_hit_intersection);
                return Lo + first_hit_intersection.vec_albedo * get_intensity(random_ray, ray_depth - 1);
            }
        }
    }
}
