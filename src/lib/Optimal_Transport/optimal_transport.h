#pragma once
#ifndef OPTIMAL_TRANSPORT_H
#define OPTIMAL_TRANSPORT_H
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <cassert>

#include "vector.h"
#include "polygon.h"
#include "lbfgs.h"
#include "arithmetic_ansi.h"
#include "particle_simulation.h"

extern std::vector<Vector> sites;
extern std::vector<double> lambdas;

lbfgsfloatval_t evaluate(
    void *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    int n,
    lbfgsfloatval_t step
);

lbfgsfloatval_t evaluate_partial(
    void *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step
);

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
);

#endif // OPTIMAL_TRANSPORT_H