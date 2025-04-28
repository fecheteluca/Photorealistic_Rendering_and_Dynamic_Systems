#include "scene.h"

namespace {
    thread_local std::default_random_engine scene_engine;
    thread_local std::uniform_real_distribution<double> scene_dist(0.0, 1.0);
    thread_local bool scene_seeded = false;
}

// Scene class definitions
Scene::Scene(
    std::vector<Geometry*> aux_l_obj, 
    Vector aux_vec_light_source, 
    Sphere* aux_sph_light_source,
    double aux_light_intensity,
    double aux_refraction_index,
    bool aux_light_source_sphere
) { 
    l_obj = aux_l_obj;
    vec_light_source = aux_vec_light_source;
    sph_light_source = aux_sph_light_source;
    light_intensity = aux_light_intensity;
    refraction_index = aux_refraction_index;
    light_source_sphere = aux_light_source_sphere;
}

void Scene::add_object(Geometry* obj_extra) {
    l_obj.push_back(obj_extra);
}

double Scene::get_random_number() {
    if (!scene_seeded) {
        #ifdef _OPENMP
            int thread_id = omp_get_thread_num();
            scene_engine.seed(15 + 1337 * thread_id);
        #else
            scene_engine.seed(15);
        #endif
            scene_seeded = true;
    }
    return scene_dist(scene_engine);
}

Intersection Scene::get_closest_hit(Ray& ray) {
    double closest_hit_distance = INT_MAX * 1.0;
    Intersection closest_hit_intersection = Intersection();
    
    for (int idx = 0; idx < l_obj.size(); idx++) {
        Intersection intersection = l_obj[idx]->intersected_by(ray);
        if (intersection.flag && intersection.distance < closest_hit_distance) {
            closest_hit_distance = intersection.distance;
            closest_hit_intersection = intersection;
            closest_hit_intersection.idx_obj = idx;
        }
    }

    return closest_hit_intersection;
}

Vector Scene::get_shadow_intensity_point_light_source(const Intersection& intersection) {
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

Vector Scene::get_shadow_intensity_sphere_light_source(const Intersection& intersection) {
    Vector vec_center = sph_light_source->get_center();
    double radius     = sph_light_source->get_radius();

    Vector vec_point  = intersection.vec_point;
    Vector vec_albedo = intersection.vec_albedo;
    Vector vec_normal = intersection.vec_normal;

    Vector vec_viz_hemishpere = (vec_point - vec_center) / (vec_point - vec_center).norm();
    Vector vec_unit_direction = random_cos(vec_viz_hemishpere);

    Vector xprime   = radius * vec_unit_direction + vec_center;
    Vector Nprime   = (xprime - vec_center) / (xprime - vec_center).norm();

    Vector omega     = (xprime - vec_point) / (xprime - vec_point).norm();
    double distance  = (xprime - vec_point).norm();

    Vector xprime_offset       = xprime + EPS * Nprime;
    Vector vec_offseted_point  = vec_point + EPS * vec_normal;

    Vector omega_offseted    = (xprime_offset - vec_offseted_point) / (xprime_offset - vec_offseted_point).norm();
    double distance_offseted = (xprime_offset - vec_offseted_point).norm();

    Ray ray(vec_offseted_point, omega_offseted);
    Intersection src_pnt_intersection = get_closest_hit(ray);

    double visibility = 0.0;
    if (distance_offseted <= src_pnt_intersection.distance) {
        visibility = 1.0;
    }

    double dot_clamped = std::max(dot(Nprime, (vec_point - vec_center) / (vec_point - vec_center).norm()), 0.0);
    double pdf = dot_clamped / (M_PI * radius * radius);

    double first_factor = (light_intensity) / (4.0 * M_PI * M_PI * radius * radius);
    Vector first_term   = first_factor * vec_albedo / M_PI;

    double second_term  = visibility
                          * std::max(dot(vec_normal, omega), 0.0)
                          * std::max(dot(Nprime, (-1) * omega), 0.0);

    double third_term   = (distance * distance) * pdf;

    return first_term * second_term / third_term;
}

Ray Scene::get_reflected_ray(Ray& ray, const Intersection& intersection) {
    Vector vec_unit_direction = ray.get_unit_direction();
    Vector vec_normal_surface = intersection.vec_normal;

    Vector vec_reflected_origin = intersection.vec_point + EPS * vec_normal_surface;

    Vector vec_reflected_unit_direction =  vec_unit_direction - 2 * dot(vec_unit_direction, vec_normal_surface) * vec_normal_surface;
    vec_reflected_unit_direction.normalize();

    Ray ray_reflected = Ray(vec_reflected_origin, vec_reflected_unit_direction);
    return ray_reflected;
}

Ray Scene::get_refracted_ray(Ray& ray, const Intersection& intersection) {
    Vector vec_unit_direction = ray.get_unit_direction();
    Vector vec_normal_surface = intersection.vec_normal;

    double cosine_unit_normal = dot(vec_unit_direction, vec_normal_surface);

    double sphere_refraction_index = l_obj[intersection.idx_obj]->get_refraction_index();
    double scene_refraction_index = refraction_index;

    double eta_1 = cosine_unit_normal <= 0 ? scene_refraction_index : sphere_refraction_index;
    double eta_2 = cosine_unit_normal <= 0 ? sphere_refraction_index : scene_refraction_index;
    double refraction_ratio = eta_1 / eta_2;

    vec_normal_surface = cosine_unit_normal <= 0 ? vec_normal_surface : (-1) * vec_normal_surface;
    cosine_unit_normal = cosine_unit_normal <= 0 ? cosine_unit_normal : (-1) * cosine_unit_normal;

    double k = 1 - pow(refraction_ratio, 2) * (1 - pow(cosine_unit_normal, 2));
            
    if (k < 0) {
        Ray ray_reflected = get_reflected_ray(ray, intersection);
        return ray_reflected;
    }
                
    Vector vec_refracted_origin = intersection.vec_point - EPS * vec_normal_surface;

    Vector vec_refracted_unit_direction_tangential = refraction_ratio * (vec_unit_direction - cosine_unit_normal * vec_normal_surface);
    Vector vec_refracted_unit_direction_normal = (-1) * vec_normal_surface * sqrt(k);
    Vector vec_refracted_unit_direction = vec_refracted_unit_direction_tangential + vec_refracted_unit_direction_normal;
    vec_refracted_unit_direction.normalize();
                    
    Ray ray_refracted = Ray(vec_refracted_origin, vec_refracted_unit_direction);
    return ray_refracted;
}

Ray Scene::get_fresnel_ray(Ray& ray, const Intersection& intersection) {
    Ray ray_reflected = get_reflected_ray(ray, intersection);
    Ray ray_refracted = get_refracted_ray(ray, intersection);

    if(ray_reflected == ray_refracted) {
        return ray_reflected;
    }

    Vector vec_unit_direction = ray.get_unit_direction();
    Vector vec_normal_surface = intersection.vec_normal;

    double cosine_unit_normal = dot(vec_unit_direction, vec_normal_surface);

    double sphere_refraction_index = l_obj[intersection.idx_obj]->get_refraction_index();
    double scene_refraction_index = refraction_index;

    double eta_1 = cosine_unit_normal <= 0 ? scene_refraction_index : sphere_refraction_index;
    double eta_2 = cosine_unit_normal <= 0 ? sphere_refraction_index : scene_refraction_index;

    double k_0 = pow((eta_1 - eta_2), 2) / pow((eta_1 + eta_2), 2);
    double R   = k_0 + (1 - k_0) * pow((1 - fabs(cosine_unit_normal)), 5);

    double u = get_random_number();
    Ray ray_returned = u < R ? ray_reflected : ray_refracted;
    return ray_returned;
}

Vector Scene::random_cos(Vector& vec) {
    double r_1 = get_random_number();
    double r_2 = get_random_number();

    double x = cos(2 * M_PI * r_1) * sqrt(1 - r_2);
    double y = sin(2 * M_PI * r_1) * sqrt(1 - r_2);;
    double z = sqrt(r_2);

    double x_vec = vec.get_x();
    double y_vec = vec.get_y();
    double z_vec = vec.get_z();

    double smallest_comp = std::min(
        fabs(x_vec),
        std::min(fabs(y_vec), fabs(z_vec))
    );    

    Vector T_1 = Vector();
    if (smallest_comp == fabs(x_vec)) {
        T_1 = Vector(0.0, z_vec, -y_vec);
    }
    else if (smallest_comp == fabs(y_vec)) {
        T_1 = Vector(-z_vec, 0.0, x_vec);
    }
    else if (smallest_comp == fabs(z_vec)) {
        T_1 = Vector(y_vec, -x_vec, 0.0);
    }
    T_1.normalize();

    Vector T_2 = cross(vec, T_1);

    return x * T_1 + y * T_2 + z * vec;
}

Ray Scene::random_ray(Intersection& intersection) {
    Vector vec_random_origin = intersection.vec_point + EPS * intersection.vec_normal;
    Vector vec_random_unit_direction = random_cos(intersection.vec_normal);
    return Ray(vec_random_origin, vec_random_unit_direction);
}

Vector Scene::get_intensity(Ray& ray , const int& ray_depth, bool last_bounce_diffuse) {
    if (ray_depth < 0) {
        return Vector(0.0, 0.0, 0.0);
    }
    else {
        Intersection first_hit_intersection = get_closest_hit(ray);
        if (first_hit_intersection.flag) {
            if (l_obj[first_hit_intersection.idx_obj]->has_mirror_surface()) {
                Ray ray_reflected = get_reflected_ray(ray, first_hit_intersection);
                return get_intensity(ray_reflected, ray_depth - 1, false);
            }
            else if (l_obj[first_hit_intersection.idx_obj]->has_transparent_surface()){
                Ray ray_returned = get_fresnel_ray(ray, first_hit_intersection);
                return get_intensity(ray_returned, ray_depth - 1, false);
            }
            else if (light_source_sphere && l_obj[first_hit_intersection.idx_obj]->is_light_source()) {
                if (last_bounce_diffuse) {
                    return Vector();
                }
                else {
                    double radius = sph_light_source->get_radius();
                    return (Vector(1.0, 1.0, 1.0) * light_intensity) / (4 * M_PI * M_PI * radius * radius);
                }
            }
            else {
                Vector Lo = Vector();
                if (!light_source_sphere) {
                    Lo = get_shadow_intensity_point_light_source(first_hit_intersection);
                }
                else {
                    Lo = get_shadow_intensity_sphere_light_source(first_hit_intersection);
                }
                Ray ray_random = random_ray(first_hit_intersection);
                return Lo + first_hit_intersection.vec_albedo * get_intensity(ray_random, ray_depth - 1, true);
            }
        }

        return Vector();
    }
}
