#ifndef SPHERE_H
#define SPHERE_H

#include <cmath>
#include "vector.h"
#include "ray.h"
#include "geometry.h"

class Sphere : public Geometry {
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

        virtual Intersection intersected_by(Ray& ray) override;

    private:
        Vector vec_center;
        double radius;
        bool invert_normals;
    };

#endif // SPHERE_H