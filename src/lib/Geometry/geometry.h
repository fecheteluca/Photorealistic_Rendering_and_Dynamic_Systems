#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "ray.h"
#include "vector.h"

struct Intersection {
    int idx_obj;
    bool flag;
    double distance;
    Vector vec_point;
    Vector vec_normal;
    Vector vec_albedo;
    

    Intersection() {
        idx_obj = -1;
        flag = false;
        distance = 0;
        vec_point = Vector();
        vec_normal = Vector();
        vec_albedo = Vector();
    }
};

class Geometry {
public:
    virtual ~Geometry() {}
    virtual Intersection intersected_by(Ray& ray) = 0;

    Vector get_color();
    bool has_mirror_surface();
    bool has_transparent_surface();
    bool is_light_source();
    double get_refraction_index();

    void set_color(const Vector& aux_vec_albedo);
    void set_mirror(const bool& aux_mirror);
    void set_transparent(const bool& aux_transparent);
    void set_light_source(const bool& aux_light_source);
    void set_refraction_index(const double& aux_refraction_index);

private:
    Vector vec_albedo;
    bool mirror;    
    bool transparent;
    bool light_source;
    double refraction_index;
};

#endif // GEOMETRY_H
