#pragma once
#ifndef PARTICLE_SIMULATION_H
#define PARTICLE_SIMULATION_H

#include <vector>
#include <ctime>
#include <random>
#include <cmath>
#include <algorithm>
#include <iostream>
#include "vector.h"
#include "polygon.h"
#include "optimal_transport.h"
#include "lbfgs.h"


// High-level fluid simulator using semi-discrete optimal transport
class ParticleSimulation {
public:
    ParticleSimulation(int N_liquid, int N_air, double dt, double eps);

    // Advance the simulation by one time step
    void Gallouet_Merigot_step();

    // Get polygons representing the current power diagram (liquid first, then air)
    std::vector<Polygon> get_current_polygons();

    int N_liquid;
    int N_air;
    double dt;
    double eps;
    double desired_fluid_volume;

    std::vector<Vector> positions_liquid;
    std::vector<Vector> velocities_liquid;
    std::vector<Vector> positions_air;
    std::vector<double> weights; // size N_liquid+1: [w_1,...,w_N_liquid, w_air]
    std::vector<double> masses;
    Vector gravity;
};

#endif // PARTICLE_SIMULATION_H
