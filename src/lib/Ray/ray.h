#ifndef RAY_H
#define RAY_H

#include "vector.h"

class Ray {
    public:
        explicit Ray(Vector aux_vec_origin = Vector(), Vector aux_vec_unit_direction = Vector());

        Vector get_origin();
        Vector get_unit_direction();

    private:
        Vector vec_origin;
        Vector vec_unit_direction;
    };

#endif // RAY_H