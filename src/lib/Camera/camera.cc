#include "camera.h"

Camera::Camera(Vector aux_vec_center, double aux_alpha, int aux_W, int aux_H) {
    vec_center = aux_vec_center;
    alpha = aux_alpha;
    W = aux_W;
    H = aux_H;
}

int Camera::get_width() {
    return W;
}

int Camera::get_height() {
    return H;
}

Ray Camera::get_ray(int i, int j) {
    // Conversion of the alpha angle to radians
    double alpha_rad = alpha * (M_PI / 180);

    // Computation of the location of pixel (i, j)
    double pixel_x = vec_center.get_x() + j + 0.5 - (W / 2);
    double pixel_y = vec_center.get_y() + (H - i - 1) + 0.5 - (H / 2);
    double pixel_z = vec_center.get_z() - (W / (2 * std::tan(alpha_rad / 2)));
    Vector vec_pixel_pos = Vector(pixel_x, pixel_y, pixel_z);

    // Computation of the direction of the ray
    Vector vec_unit_direction = vec_pixel_pos - vec_center;
    vec_unit_direction.normalize();

    return Ray(vec_center, vec_unit_direction);
}