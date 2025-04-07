#include "ray.h"

// Ray class definitions
Ray::Ray(
    Vector aux_vec_origin, 
    Vector aux_vec_unit_direction
) {
    vec_origin = aux_vec_origin;
    vec_unit_direction = aux_vec_unit_direction;
}

Vector Ray::get_origin() {
    return vec_origin;
}

Vector Ray::get_unit_direction() {
    return vec_unit_direction;
}