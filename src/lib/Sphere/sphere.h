#ifndef SPHERE_H
#define SPHERE_H

#include <cmath>
#include "vector.h"
#include "ray.h"

struct Intersection {
    bool flag;
    double distance;
    Vector vec_point;
    Vector vec_normal;
    Vector vec_albedo;

    Intersection() {
        flag = false;
        distance = 0;
        vec_point = Vector();
        vec_normal = Vector();
        vec_albedo = Vector();
    }
};

class Sphere {
    public:
        explicit Sphere(Vector aux_vec_center = Vector(), double aux_radius = 0, Vector aux_vec_albedo = Vector());

        Vector get_center();
        double get_radius();
        Vector get_color();

        Intersection intersected_by(Ray ray);

    private:
        Vector vec_center;
        double radius;
        Vector vec_albedo;
    };

#endif // SPHERE_H