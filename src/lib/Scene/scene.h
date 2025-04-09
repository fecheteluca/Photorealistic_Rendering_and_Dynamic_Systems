#ifndef SCENE_H
#define SCENE_H

#include <iostream>
#include <cstdlib>
#include <random>
#include <climits>
#include <vector>
#include <cmath>
#include <omp.h>

#include "vector.h"
#include "ray.h"
#include "sphere.h"

#define EPS 0.001

class Scene {
    public:
        explicit Scene(
            std::vector<Sphere> aux_l_sph = std::vector<Sphere>(), 
            Vector aux_vec_light_source = Vector(), 
            double aux_light_intensity = 0,
            double aux_refraction_index = 1.0
        );

        void add_object(const Sphere& sph_extra);

        Intersection get_closest_hit(Ray& ray);
        
        Vector get_shadow_intensity(const Intersection& intersection);

        Ray random_cos(const Intersection& intersection);
        Vector get_intensity(Ray& ray, const int& ray_depth);
        
    private:
        std::vector<Sphere> l_sph;
        Vector vec_light_source;
        double light_intensity;
        double refraction_index;
    };

#endif // SCENE_H