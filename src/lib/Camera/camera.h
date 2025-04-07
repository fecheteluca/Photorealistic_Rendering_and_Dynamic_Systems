#ifndef CAMERA_H
#define CAMERA_H

#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <cmath>
#include "vector.h"
#include "ray.h"
#include "sphere.h"
#include "scene.h"

class Camera {
    public:
        explicit Camera(
            Vector aux_vec_center = Vector(), 
            double aux_alpha = 0, 
            int aux_W = 0, 
            int aux_H = 0
        );
        
        int get_width();
        int get_height();

        Ray get_ray(int i, int j);

    private:
        Vector vec_center; // The origin of the camera
        double alpha; // The visual angle covering the W pixels in width
        int W; // The width of the rendered image
        int H; // The height of the rendered image
    };

#endif // CAMERA_H
