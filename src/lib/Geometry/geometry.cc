#include "geometry.h"


Vector Geometry::get_color() {
    return vec_albedo;
}

bool Geometry::has_mirror_surface() {
    return mirror;
}

bool Geometry::has_transparent_surface() {
    return transparent;
}

bool Geometry::is_light_source() {
    return light_source;
}

double Geometry::get_refraction_index() {
    return refraction_index;
}

void Geometry::set_color(const Vector& aux_vec_albedo) {
    vec_albedo = aux_vec_albedo;
}

void Geometry::set_mirror(const bool& aux_mirror) {
    mirror = aux_mirror;
}

void Geometry::set_transparent(const bool& aux_transparent) {
    transparent = aux_transparent;
}

void Geometry::set_light_source(const bool& aux_light_source) {
    light_source = aux_light_source;
}

void Geometry::set_refraction_index(const double& aux_refraction_index) {
    refraction_index = aux_refraction_index;
}
