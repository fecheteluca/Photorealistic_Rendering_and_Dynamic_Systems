#ifndef SCENE_H
#define SCENE_H

#include <vector>
#include <climits>

#include "vector.h"
#include "ray.h"
#include "sphere.h"

class Scene {
    public:
        explicit Scene(std::vector<Sphere> aux_l_sph = std::vector<Sphere>());

        void add_object(Sphere sph_extra);

        Intersection get_closest_hit(Ray ray);
        
    private:
        std::vector<Sphere> l_sph;
    };

    Scene get_standard_scene();

#endif // SCENE_H