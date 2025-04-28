#include "sphere.h"

// Sphere class definitions
Sphere::Sphere(
    Vector aux_vec_center, 
    double aux_radius, 
    Vector aux_vec_albedo, 
    bool aux_mirror,
    bool aux_transparent,
    bool aux_light_source,
    double aux_refraction_index,
    bool aux_invert_normals
) {
    vec_center = aux_vec_center;
    radius = aux_radius;
    invert_normals = aux_invert_normals;

    this->set_color(aux_vec_albedo);
    this->set_mirror(aux_mirror);
    this->set_transparent(aux_transparent);
    this->set_light_source(aux_light_source);
    this->set_refraction_index(aux_refraction_index);
}

Vector Sphere::get_center() {
    return vec_center;
}

double Sphere::get_radius() {
    return radius;
}

Intersection Sphere::intersected_by(Ray& ray) {
    Intersection intersection = Intersection();
    
    Vector vec_origin = ray.get_origin();
    Vector vec_unit_direction = ray.get_unit_direction();

    Vector vec_org_cnt_diff         = vec_origin - vec_center;
    double discriminant_first_elem  = pow(dot(vec_unit_direction, vec_org_cnt_diff), 2);
    double discriminant_second_elem = pow(vec_org_cnt_diff.norm(), 2) - pow(radius, 2);
    double discriminant             = discriminant_first_elem - discriminant_second_elem;

    if (discriminant >= 0) {
        double first_intersection  = dot(vec_unit_direction, vec_center - vec_origin) - sqrt(discriminant);
        double second_intersection = dot(vec_unit_direction, vec_center - vec_origin) + sqrt(discriminant);
        
        if (second_intersection >= 0 && first_intersection >= 0) {
            intersection.flag = true;
            intersection.distance = first_intersection;
            intersection.vec_point = vec_origin + first_intersection * vec_unit_direction;
            intersection.vec_normal = !invert_normals ? intersection.vec_point - vec_center : vec_center - intersection.vec_point;
            intersection.vec_normal.normalize();
            intersection.vec_albedo = this->get_color();
        }
        else if (second_intersection >= 0 && first_intersection < 0) {
            intersection.flag = true;
            intersection.distance = second_intersection;
            intersection.vec_point = vec_origin + second_intersection * vec_unit_direction;
            intersection.vec_normal = !invert_normals ? intersection.vec_point - vec_center : vec_center - intersection.vec_point;
            intersection.vec_normal.normalize();
            intersection.vec_albedo = this->get_color();
        }
    }

    return intersection;
}