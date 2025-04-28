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
#include "geometry.h"

#define EPS 0.001

class Scene {
    public:
        explicit Scene(
            std::vector<Geometry*> aux_l_sph = std::vector<Geometry*>(), 
            Vector aux_vec_light_source = Vector(), 
            Sphere* sph_light_source = new Sphere(),
            double aux_light_intensity = 0,
            double aux_refraction_index = 1.0,
            bool aux_light_source_sphere = false
        );

        void add_object(Geometry* obj_extra);

        double get_random_number();

        Intersection get_closest_hit(Ray& ray);
        
        Vector get_shadow_intensity_point_light_source(const Intersection& intersection);
        Vector get_shadow_intensity_sphere_light_source(const Intersection& intersection);

        Ray get_reflected_ray(Ray& ray, const Intersection& intersection);
        Ray get_refracted_ray(Ray& ray, const Intersection& intersection);
        Ray get_fresnel_ray(Ray& ray, const Intersection& intersection);

        Vector random_cos(Vector& vec);
        Ray random_ray(Intersection& intersection);

        Vector get_intensity(Ray& ray, const int& ray_depth, bool last_bounce_diffuse);
        
    private:
        std::vector<Geometry*> l_obj;
        Vector vec_light_source;
        Sphere* sph_light_source;
        double light_intensity;
        double refraction_index;
        bool light_source_sphere;
    };

#endif // SCENE_H