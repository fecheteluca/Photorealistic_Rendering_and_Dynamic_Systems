#ifndef SCENE_H
#define SCENE_H

#include <vector>
#include <climits>

#include "vector.h"
#include "ray.h"
#include "sphere.h"

class Scene {
    public:
        explicit Scene(std::vector<Sphere> aux_l_sph = std::vector<Sphere>(), Vector aux_vec_light_source = Vector(), double aux_light_intensity = 0);

        void add_object(Sphere sph_extra);

        Intersection get_closest_hit(Ray ray);

        Vector get_shadow_intensity(Intersection intersection);
        
    private:
        std::vector<Sphere> l_sph;
        Vector vec_light_source;
        double light_intensity;
    };

    Scene get_standard_scene();

#endif // SCENE_H