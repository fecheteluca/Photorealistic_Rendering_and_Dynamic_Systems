#ifndef SCENE_H
#define SCENE_H

#include <iostream>
#include <climits>
#include <vector>

#include "vector.h"
#include "ray.h"
#include "sphere.h"

class Scene {
    public:
        explicit Scene(
            std::vector<Sphere> aux_l_sph = std::vector<Sphere>(), 
            Vector aux_vec_light_source = Vector(), 
            double aux_light_intensity = 0,
            double aux_refraction_index = 1.0
        );

        void add_object(Sphere sph_extra);

        Intersection get_closest_hit(Ray ray);

        Vector get_shadow_intensity(Intersection intersection);

        Vector get_reflected_intensity(Ray ray , int ray_depth);

        Vector get_refracted_intensity(Ray ray, int ray_depth);
        
    private:
        std::vector<Sphere> l_sph;
        Vector vec_light_source;
        double light_intensity;
        double refraction_index;
    };

    Scene get_standard_scene();
    Scene get_intermediate_scene_reflection();
    Scene get_intermediate_scene_refraction_reflection();

#endif // SCENE_H