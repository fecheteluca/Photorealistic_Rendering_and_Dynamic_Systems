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
    double aux_spread,
    double aux_focal_dist,
    double aux_aperture
) { 
    vec_center = aux_vec_center;
    alpha      = aux_alpha;
    W          = aux_W;
    H          = aux_H;
    stdev      = aux_stdev;
    spread     = aux_spread;
    focal_dist = aux_focal_dist;
    aperture   = aux_aperture;
}

int Camera::get_width() {
    return W;
}

int Camera::get_height() {
    return H;
}

double Camera::get_random_number() {
    if (!camera_seeded) {
        #ifdef _OPENMP
            int thread_id = omp_get_thread_num();
            camera_engine.seed(15 + 1337 * thread_id);
        #else
            camera_engine.seed(15);
        #endif
            camera_seeded = true;
    }

    return camera_dist(camera_engine);
}

void Camera::boxMuller(double stdev , double &x , double &y) {
    double r1 = get_random_number();
    double r2 = get_random_number();
    x = sqrt((-2) * log (r1)) * cos (2 * M_PI * r2) * stdev;
    y = sqrt((-2) * log (r1)) * sin (2 * M_PI * r2) * stdev;
}

Ray Camera::get_ray(const int& i, const int& j) {
    double alpha_rad = alpha * (M_PI / 180.0);

    double dx, dy;
    boxMuller(stdev, dx, dy); 
    dx *= spread;
    dy *= spread;

    // Computation of the location of pixel (i, j)
    double pixel_x = vec_center.get_x() + (j + dx) + 0.5 - (W / 2.0);
    double pixel_y = vec_center.get_y() + (H - (i + dy) - 1) + 0.5 - (H / 2.0);
    double pixel_z = vec_center.get_z() - (W / (2 * std::tan(alpha_rad / 2)));
    Vector vec_pixel_pos = Vector(pixel_x, pixel_y, pixel_z);

    // Computation of the direction of the ray
    Vector vec_unit_direction = vec_pixel_pos - vec_center;
    vec_unit_direction.normalize();

    if (aperture <= 1e-8) {
        return Ray(vec_center, vec_unit_direction);
    }

    double t_focus = focal_dist / std::fabs(vec_unit_direction.get_z());
    Vector vec_focus_point = vec_center + (vec_unit_direction * t_focus);

    double r     = std::sqrt(get_random_number());
    double theta = 2.0 * M_PI * get_random_number();
    double lens_r = aperture * r;

    double lens_x = lens_r * std::cos(theta);
    double lens_y = lens_r * std::sin(theta);

    Vector vec_dof_origin = Vector(
        vec_center.get_x() + lens_x,
        vec_center.get_y() + lens_y,
        vec_center.get_z()
    );

    Vector vec_dof_dir = vec_focus_point - vec_dof_origin;
    vec_dof_dir.normalize();

    return Ray(vec_dof_origin, vec_dof_dir);
}