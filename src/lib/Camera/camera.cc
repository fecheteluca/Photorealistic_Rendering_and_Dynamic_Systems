#include "camera.h"

namespace {
    thread_local std::default_random_engine camera_engine;
    thread_local std::uniform_real_distribution<double> camera_dist(0.0, 1.0);
    thread_local bool camera_seeded = false;
}

Camera::Camera(
    Vector aux_vec_center, 
    double aux_alpha, 
    int aux_W, 
    int aux_H,
    double aux_stdev,
    double aux_spread
) { 
    vec_center = aux_vec_center;
    alpha = aux_alpha;
    W = aux_W;
    H = aux_H;
    stdev = aux_stdev;
    spread = aux_spread;
}

int Camera::get_width() {
    return W;
}

int Camera::get_height() {
    return H;
}

void Camera::boxMuller(double stdev , double &x , double &y) {
    if (!camera_seeded) {
        #ifdef _OPENMP
            int thread_id = omp_get_thread_num();
            camera_engine.seed(15 + 1337 * thread_id);
        #else
            camera_engine.seed(15);
        #endif
            camera_seeded = true;
    }
    
    double r1 = camera_dist(camera_engine);
    double r2 = camera_dist(camera_engine);
    x = sqrt((-2) * log (r1)) * cos (2 * M_PI * r2) * stdev;
    y = sqrt((-2) * log (r1)) * sin (2 * M_PI * r2) * stdev;
}

Ray Camera::get_ray(const int& i, const int& j) {
    // Conversion of the alpha angle to radians
    double alpha_rad = alpha * (M_PI / 180);

    double dx, dy;
    boxMuller(stdev, dx, dy); 
    dx *= spread;
    dy *= spread;

    // Computation of the location of pixel (i, j)
    double pixel_x = vec_center.get_x() + (j + dy) + 0.5 - (W / 2.0);
    double pixel_y = vec_center.get_y() + (H - (i + dx) - 1) + 0.5 - (H / 2.0);
    double pixel_z = vec_center.get_z() - (W / (2 * std::tan(alpha_rad / 2)));
    Vector vec_pixel_pos = Vector(pixel_x, pixel_y, pixel_z);

    // Computation of the direction of the ray
    Vector vec_unit_direction = vec_pixel_pos - vec_center;
    vec_unit_direction.normalize();

    return Ray(vec_center, vec_unit_direction);
}