#include "particle_simulation.h"

ParticleSimulation::ParticleSimulation(int N_liquid_, int N_air_, double dt_, double eps_)
    : N_liquid(N_liquid_), N_air(N_air_), dt(dt_), eps(eps_), gravity(0, -2.0)
{
    Vector center = Vector(0.5, 0.7);
    double radius = 0.18;
    positions_liquid.reserve(N_liquid);
    for (int i = 0; i < N_liquid; ++i) {
        double r = radius * std::sqrt(std::rand() / (double)RAND_MAX);
        double theta = 2.0 * M_PI * (std::rand() / (double)RAND_MAX);
        positions_liquid.push_back(Vector(center[0] + r * std::cos(theta), center[1] + r * std::sin(theta)));
    }

    positions_air.reserve(N_air);
    for (int i = 0; i < N_air; ++i) {
        positions_air.push_back(Vector(std::rand() / (double)RAND_MAX, std::rand() / (double)RAND_MAX));
    }

    velocities_liquid.assign(N_liquid, Vector(0,0));

    const int steps = 20;
    for (int iter = 0; iter < steps; iter++) {
        std::vector<Polygon> cells_liquid = voronoi_ple(positions_liquid);
        for (int i = 0; i < N_liquid; ++i) {
            positions_liquid[i] = compute_moments(clip_by_disk(cells_liquid[i], center, radius)).centroid;
        }
        
        std::vector<Polygon> cells_air = voronoi_ple(positions_air);
        for (int j = 0; j < N_air; ++j) {
            positions_air[j] = compute_moments(cells_air[j]).centroid;
        }
    }

    desired_fluid_volume = M_PI * radius * radius;

    weights.resize(N_liquid + 1, 0.5);
    masses.resize(N_liquid, 2.0);
}

void ParticleSimulation::Gallouet_Merigot_step() {
    // Optimizing the weights for each Gallouet Merigot step
    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);
    param.max_iterations    = 500;
    param.epsilon           = 1e-6;
    param.past              = 3;
    param.delta             = 1e-6;
    param.max_linesearch    = 100;
    param.max_linesearch    = 200;
    param.max_step          = 1e8;  
    param.min_step          = 1e-20;
    param.linesearch        = LBFGS_LINESEARCH_BACKTRACKING;

    std::vector<lbfgsfloatval_t> x(weights.begin(), weights.end());
    lbfgsfloatval_t final_F;
    int ret = lbfgs(
        N_liquid + 1, 
        x.data(), 
        &final_F,
        evaluate_partial,
        lbfgs_progress,
        this,
        &param
    );

    weights.assign(x.begin(), x.end());

    std::vector<Vector> all_sites;
    for (int i = 0; i < N_liquid; i++) all_sites.push_back(positions_liquid[i]);
    for (int i = 0; i < N_air; i++) all_sites.push_back(positions_air[i]);

    std::vector<double> all_weights;
    for (int i = 0; i < N_liquid; i++) all_weights.push_back(weights[i]);
    for (int i = 0; i < N_air; i++) all_weights.push_back(weights[N_liquid]);

    auto cells = weighted_voronoi_ple(all_sites, all_weights);

    for (int i = 0; i < N_liquid; i++) {
        double rad = std::sqrt(std::max(0.0, weights[i] - weights[N_liquid]));
        Vector centroid = compute_moments(clip_by_disk(cells[i], positions_liquid[i], rad)).centroid;

        Vector spring_force = (1.0 / (eps * eps)) * (centroid - positions_liquid[i]);
        Vector total_force = spring_force + masses[i] * gravity;
        
        velocities_liquid[i] = velocities_liquid[i] + (dt / masses[i]) * total_force;
        positions_liquid[i] = positions_liquid[i] +  dt * velocities_liquid[i];

        // If the liquid particles go outside the domain (unit square), they bounce back
        if (positions_liquid[i][0] < 0) {
            positions_liquid[i][0] = 0;
            velocities_liquid[i][0] *= -0.5;
        }
        if (positions_liquid[i][0] > 1) {
            positions_liquid[i][0] = 1;
            velocities_liquid[i][0] *= -0.5;
        }
        if (positions_liquid[i][1] < 0) {
            positions_liquid[i][1] = 0;
            velocities_liquid[i][1] *= -0.5;
        }
        if (positions_liquid[i][1] > 1) {
            positions_liquid[i][1] = 1;
            velocities_liquid[i][1] *= -0.5;
        }
    }
}

std::vector<Polygon> ParticleSimulation::get_current_polygons() {
    std::vector<Vector> all_sites;
    for (int i = 0; i < N_liquid; i++) all_sites.push_back(positions_liquid[i]);
    for (int i = 0; i < N_air; i++) all_sites.push_back(positions_air[i]);

    std::vector<double> all_weights;
    for (int i = 0; i < N_liquid; i++) all_weights.push_back(weights[i]);
    for (int i = 0; i < N_air; i++) all_weights.push_back(weights[N_liquid]);

    std::vector<Polygon> cells = weighted_voronoi_ple(all_sites, all_weights);

    std::vector<Polygon> result;
    for (int i = 0; i < N_liquid; ++i) {
        double rad = std::sqrt(std::max(0.0, weights[i] - weights[N_liquid]));
        result.push_back(clip_by_disk(cells[i], positions_liquid[i], rad));
    }
    for (int i = N_liquid; i < N_liquid + N_air; ++i) {
        result.push_back(cells[i]);
    }
    
    return result;
}
