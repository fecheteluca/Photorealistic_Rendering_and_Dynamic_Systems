#include "optimal_transport.h"

// the global problem data for evaluate()
std::vector<Vector> sites;
std::vector<double> lambdas;

// Semi-Discrete Optimal Transport
lbfgsfloatval_t evaluate(
    void * instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step
) {
    std::vector<double> weights(n);
    for (int i = 0; i < n; ++i) weights[i] = static_cast<double>(x[i]);

    auto cells = weighted_voronoi_ple(sites, weights);

    double gW_value = 0.0;
    for (int i = 0; i < n; ++i) {
        Polygon &polygon = cells[i];

        double area = compute_moments(polygon).signed_area;

        double integral_over_polygon = 0.0;
        Vector &fixed_vertex = polygon.vertices[0]; // We fix a vertex and we build triangles to compute the integral term
        for (int k = 1; k < (int)polygon.vertices.size() - 1; k++) {
            Vector &vertex_k = polygon.vertices[k];
            Vector &vertex_k_plus_1 = polygon.vertices[k + 1];

            // The value of the integral over a triangle of P -> \|P - P_i\|^2
            double kth_triangle_area = compute_moments(Polygon({fixed_vertex, vertex_k, vertex_k_plus_1})).signed_area;
            integral_over_polygon += (kth_triangle_area / 6.0) * (
                     dot(fixed_vertex - sites[i], fixed_vertex - sites[i])
                   + dot(fixed_vertex - sites[i], vertex_k - sites[i])
                   + dot(fixed_vertex - sites[i], vertex_k_plus_1 - sites[i])
                   + dot(vertex_k - sites[i], vertex_k - sites[i])
                   + dot(vertex_k - sites[i], vertex_k_plus_1 - sites[i])
                   + dot(vertex_k_plus_1 - sites[i], vertex_k_plus_1 - sites[i]));
        }
        // Because the function f is constant, the integral over T of f equals the area of T times f
        gW_value += integral_over_polygon - weights[i] * area + lambdas[i] * weights[i];
        g[i] = area - lambdas[i]; // The i-th value of the gradient of g
    }
    return static_cast<lbfgsfloatval_t>(-gW_value);
}

lbfgsfloatval_t evaluate_partial(
    void * instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step
) {
    ParticleSimulation* sim = static_cast<ParticleSimulation*>(instance);
    const int N_liquid = sim->N_liquid;
    const int N_air = sim->N_air;

    std::vector<double> w_liquid;
    for (int i = 0; i < N_liquid; i++) 
        w_liquid.push_back(static_cast<double>(x[i]));
    const double w_air = x[N_liquid];

    std::vector<double> all_weights;
    for (int i = 0; i < N_liquid; i++) all_weights.push_back(w_liquid[i]);
    for (int i = 0; i < N_air; i++) all_weights.push_back(w_air);

    std::vector<Vector> all_sites;
    for (int i = 0; i < N_liquid; i++) all_sites.push_back(sim->positions_liquid[i]);
    for (int i = 0; i < N_air; i++) all_sites.push_back(sim->positions_air[i]);

    std::vector<double> all_lambdas;
    for (int i = 0; i < N_liquid; i++) all_lambdas.push_back(sim->desired_fluid_volume / N_liquid);
    for (int i = 0; i < N_air; i++) all_lambdas.push_back((1 - sim->desired_fluid_volume) / N_air);

    std::vector<Polygon> cells = weighted_voronoi_ple(all_sites, all_weights);

    double desired_fluid_volume = sim->desired_fluid_volume;
    double desired_air_volume = 1.0 - desired_fluid_volume;
    double total_fluid_area = 0.0;

    double gW_value = 0.0;
    for (int i = 0; i < N_liquid; i++) {
        Polygon &polygon = cells[i];

        double area = compute_moments(polygon).signed_area;
        total_fluid_area += compute_moments(polygon).area;

        double integral_over_polygon = 0.0;
        Vector &fixed_vertex = polygon.vertices[0]; 
        for (int k = 1; k < (int)polygon.vertices.size() - 1; k++) {
            Vector &vertex_k = polygon.vertices[k];
            Vector &vertex_k_plus_1 = polygon.vertices[k + 1];

            double kth_triangle_area = compute_moments(Polygon({fixed_vertex, vertex_k, vertex_k_plus_1})).signed_area;
            integral_over_polygon += (kth_triangle_area / 6.0) * (
                     dot(fixed_vertex - all_sites[i], fixed_vertex - all_sites[i])
                   + dot(fixed_vertex - all_sites[i], vertex_k - all_sites[i])
                   + dot(fixed_vertex - all_sites[i], vertex_k_plus_1 - all_sites[i])
                   + dot(vertex_k - all_sites[i], vertex_k - all_sites[i])
                   + dot(vertex_k - all_sites[i], vertex_k_plus_1 - all_sites[i])
                   + dot(vertex_k_plus_1 - all_sites[i], vertex_k_plus_1 - all_sites[i]));
        }
        gW_value += integral_over_polygon - all_weights[i] * area + all_lambdas[i] * all_weights[i];
        g[i] = all_lambdas[i] - area; 
    }

    double estimated_air_volume = 1.0 - total_fluid_area;
    gW_value += w_air * (desired_air_volume - estimated_air_volume);
    g[N_liquid] = desired_air_volume / N_air  - estimated_air_volume; 

    return static_cast<lbfgsfloatval_t>(gW_value);
}

int lbfgs_progress(
    void * instance,
    const lbfgsfloatval_t *x,
    const lbfgsfloatval_t *g,
    const lbfgsfloatval_t fx,
    const lbfgsfloatval_t xnorm,
    const lbfgsfloatval_t gnorm,
    const lbfgsfloatval_t step,
    int n,
    int k,
    int ls
) {
    std::cout
      << "Iter " << k
      << " | F="   << -fx
      << " | ||g||=" << gnorm
      << " | step="  << step
      << std::endl;
    return 0;
}
