#include "camera.h"

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
    double r1 = uniform ( engine ) ;
    double r2 = uniform ( engine ) ;
    x = sqrt((-2) * log ( r1 ) ) *cos ( 2 * M_PI*r2 ) *stdev ;
    y = sqrt((-2) * log ( r1 ) ) *sin ( 2 * M_PI*r2 ) *stdev ;
}

Ray Camera::get_ray(const int& i, const int& j) {
    // Conversion of the alpha angle to radians
    double alpha_rad = alpha * (M_PI / 180);

    double dx, dy;
    boxMuller(stdev, dx, dy); 
    dx *= spread;
    dy *= spread;

    // Computation of the location of pixel (i, j)
    double pixel_x = vec_center.get_x() + (j + 0.5 + dx) - (W / 2.0);
    double pixel_y = vec_center.get_y() + ((H - i - 1) + 0.5 + dy) - (H / 2.0);
    double pixel_z = vec_center.get_z() - (W / (2 * std::tan(alpha_rad / 2)));
    Vector vec_pixel_pos = Vector(pixel_x, pixel_y, pixel_z);

    // Computation of the direction of the ray
    Vector vec_unit_direction = vec_pixel_pos - vec_center;
    vec_unit_direction.normalize();

    return Ray(vec_center, vec_unit_direction);
}