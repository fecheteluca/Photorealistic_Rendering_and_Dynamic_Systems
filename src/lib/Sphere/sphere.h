#ifndef SPHERE_H
#define SPHERE_H

#include <cmath>
#include "vector.h"
#include "ray.h"

struct Intersection {
    int idx_sph;
    bool flag;
    double distance;
    Vector vec_point;
    Vector vec_normal;
    Vector vec_albedo;
    

    Intersection() {
        idx_sph = -1;
        flag = false;
        distance = 0;
        vec_point = Vector();
        vec_normal = Vector();
        vec_albedo = Vector();
    }
};

class Sphere {
    public:
        explicit Sphere(
            Vector aux_vec_center = Vector(), 
            double aux_radius = 0, 
            Vector aux_vec_albedo = Vector(), 
            bool aux_mirror = false, 
            bool aux_transparent = false,
            bool aux_light_source = false,
            double aux_refraction_index = 1.0,
            bool aux_invert_normals = false
        );

        Vector get_center();
        double get_radius();
        Vector get_color();
        bool has_mirror_surface();
        bool has_transparent_surface();
        bool is_light_source();
        double get_refraction_index();

        Intersection intersected_by(Ray& ray);

    private:
        Vector vec_center;
        double radius;
        Vector vec_albedo;
        bool mirror;    
        bool transparent;
        bool light_source;
        double refraction_index;
        bool invert_normals;
    };

#endif // SPHERE_H